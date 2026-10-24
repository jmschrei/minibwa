"""Area A: byte pins of `minibwa map` output on inputs that stress seeding.

The other area A tests check intervals and anchors against brute force; these
pin what the whole mapping path makes of them, so that a change in seeding
that alters any record (a sampled anchor, a sub-SMEM, a batch boundary) shows
up as a diff. The inputs are simulated from the standard synthetic genome
(conftest.make_synth: repeats, low-complexity tracts, N runs, short contigs):
short reads from anywhere in it, reads with N, IUPAC codes and lowercase
bases, reads across contig junctions and a few long reads. Pins are created
on base minibwa with MINIBWA_TEST_UPDATE_PINS=1.
"""

import pytest

import minibwa_helpers as H
import helpers_a as A


def _reads(genome):
	rng = H.Rng(9001)
	out = list(H.simulate_reads(genome, 360, length=(60, 250), seed=51, sub_rate=0.01, indel_rate=0.001,
		n_rate=0.002, name_prefix="sr"))
	out += H.simulate_reads(genome, 40, length=150, seed=52, sub_rate=0.0, indel_rate=0.0, contigs=None,
		name_prefix="ex")
	# reads with IUPAC codes, '.', and lowercase stretches
	for i, r in enumerate(H.simulate_reads(genome, 30, length=150, seed=53, name_prefix="amb")):
		s = list(r.seq)
		for _ in range(1 + i % 4):
			s[rng.below(len(s))] = rng.choice("RYKMN.n")
		if i % 3 == 0:
			a = rng.below(100)
			s[a:a + 40] = [c.lower() for c in s[a:a + 40]]
		out.append((r.name, "".join(s), r.qual))
	# across the junctions between consecutive contigs, both orientations
	for i in range(len(genome.names) - 1):
		a, b = genome[genome.names[i]], genome[genome.names[i + 1]]
		for k, (la, lb) in enumerate(((70, 80), (40, 110), (120, 30))):
			s = a[-la:] + b[:lb]
			out.append(("junc%d_%d" % (i, k), s if k % 2 == 0 else H.revcomp(s), "I" * len(s)))
	# low-complexity and repeat copies
	for j, f in enumerate(genome.features):
		if f.kind in ("lowcomp", "repeat"):
			st = max(0, f.start - 40)
			s = genome[f.contig][st:f.end + 40][:300]
			out.append(("feat%d_%s" % (j, f.kind), s if j % 2 else H.revcomp(s), "I" * len(s)))
	out += H.simulate_long_reads(genome, 6, length=(1500, 3000), seed=54)
	return out


@pytest.fixture(scope="module")
def pin_reads(synth, tmp_path_factory):
	return H.write_fastq(str(tmp_path_factory.mktemp("pins") / "reads.fq"), _reads(synth.genome))


@pytest.mark.pinned
@pytest.mark.parametrize("name,args", [
	("a_synth_default.sam.gz", []),
	("a_synth_k15_c5.sam.gz", ["-k", "15", "-c", "5"]),
	("a_synth_k31_c1000.sam.gz", ["-k", "31", "-c", "1000", "-b", "cs"]),
	("a_synth_k12_c2_paf.gz", ["-k", "12", "-c", "2", "-f"]),
	("a_synth_sr_k21.sam.gz", ["-x", "sr", "-k", "21"]),
])
def test_pinned_map_output(synth, pin_reads, threads, name, args):
	"""Pins the current output of minibwa map on 448 short and 6 long single-end reads (see
	the module docstring) under several seeding options; every read has one primary record
	consistent with the genome. The output does not depend on -t."""
	r = H.run_map(synth.index, pin_reads, args)
	text = H.strip_pg(r.stdout)
	if "-f" not in args:
		sam = H.parse_sam(text)
		assert len(sam.primaries()) == len(H.read_fastx(pin_reads))
		assert H.sam_invariant_errors(sam, synth.genome) == []
	H.check_pinned(name, text)
