"""Area C: byte pins of `minibwa map` paired-end and Hi-C output on fixed synthetic
input, for the code paths a speed optimization of pairing is likely to touch:
insert-size inference, pairing among repeat copies, mate rescue, pair MAPQ,
-5 and --hic, and libraries in the FF, RR and RF orientations as well as FR.

The genome has diverged repeat families (so pairing decides between copies
and the pair MAPQ depends on the insert-size statistics) and the reads
include chimeric read 1s (so -5 changes primaries), mates without seeds
(rescue) and random mates. Each test pins the current output of base minibwa
(74bdd48); regenerate with MINIBWA_TEST_UPDATE_PINS=1 only on base minibwa.
"""

import pytest

import minibwa_helpers as H
import helpers_c as C


def pin_genome():
	return H.make_genome(seed=95, contig_lengths=(80000, 50000), repeat_families=8, repeat_len=(200, 1200),
		repeat_copies=4, repeat_divergence=0.02, low_complexity=3, n_runs=2, short_contigs=1)


def _pe_reads(g):
	"""400 FR pairs from anywhere, 25 pairs with a chimeric read 1 (a 40-70 bp 5' part from
	elsewhere), 20 pairs whose read 2 has no seed, and 5 pairs with a random read 2."""
	pairs = H.simulate_pairs(g, 400, seed=91, insert_mean=320, insert_sd=35, sub_rate=0.01, indel_rate=0.001,
		qual=None, name_prefix="p")
	recs1 = [(p.name, p.r1.seq) for p in pairs]
	recs2 = [(p.name, p.r2.seq) for p in pairs]
	rng = H.Rng(96)
	chim = H.simulate_pairs(g, 25, seed=97, insert_mean=320, insert_sd=35, sub_rate=0.005, indel_rate=0.0,
		qual=None, name_prefix="c")
	other = H.simulate_reads(g, 25, length=(40, 70), seed=98, sub_rate=0.005, indel_rate=0.0, qual=None)
	for p, o in zip(chim, other):
		recs1.append((p.name, o.seq + p.r1.seq[len(o.seq):]))
		recs2.append((p.name, p.r2.seq))
	resc = H.simulate_pairs(g, 20, seed=92, insert_mean=320, insert_sd=35, sub_rate=0.0, indel_rate=0.0,
		qual=None, name_prefix="s", unique=True)
	for p in resc:
		recs1.append((p.name, p.r1.seq))
		recs2.append((p.name, C.substitute(p.r2.seq, C.unseedable_subs(len(p.r2.seq)))))
	junk = C.PlantedGenome(93, [("x", 750)])["x"]
	for i in range(5):
		recs1.append(("j%d" % i, pairs[i].r1.seq))
		recs2.append(("j%d" % i, junk[i * 150:(i + 1) * 150]))
	return recs1, recs2


@pytest.fixture(scope="module")
def pin_input(index_of, tmp_path_factory):
	g = pin_genome()
	idx = index_of(g)
	d = tmp_path_factory.mktemp("cpins")
	recs1, recs2 = _pe_reads(g)
	f1 = H.write_fasta(str(d / "r1.fa"), recs1)
	f2 = H.write_fasta(str(d / "r2.fa"), recs2)
	hic = H.simulate_hic_pairs(g, 150, seed=94, chimeric_frac=0.4, unique=False, qual=None)
	h1, h2 = H.write_pairs(str(d / "h1.fa"), str(d / "h2.fa"), hic)
	return g, idx, f1, f2, h1, h2


PE_RUNS = [
	("c_pe_default.sam.gz", []),
	("c_pe_primary5.sam.gz", ["-5"]),
	("c_pe_fixed.sam.gz", ["-I", "300,60", "--rescue=3", "-U", "25"]),
]


@pytest.mark.pinned
@pytest.mark.parametrize("name,args", PE_RUNS, ids=["default", "primary5", "fixed"])
def test_pe_output_pinned(pin_input, name, args):
	"""Pins the current output: paired-end SAM with inferred (default, -5) or fixed insert-size
	statistics, pairing among repeat copies, chimeric read 1s, mate rescue and unmapped mates."""
	g, idx, f1, f2, h1, h2 = pin_input
	res = H.run_map(idx, [f1, f2], args)
	sam = res.sam()
	assert H.pair_invariant_errors(sam) == []
	assert sum(r.is_proper_pair for r in sam.records if r.is_primary) > 600
	H.check_pinned(name, H.strip_pg(res.stdout))


@pytest.mark.pinned
def test_pe_paf_pinned(pin_input):
	"""Pins the current output: paired-end PAF (names with /1 and /2)."""
	g, idx, f1, f2, h1, h2 = pin_input
	res = H.run_map(idx, [f1, f2], ["-f"])
	assert all(p.qname[-2:] in ("/1", "/2") for p in res.paf())
	H.check_pinned("c_pe.paf.gz", res.stdout)


@pytest.mark.pinned
def test_hic_output_pinned(pin_input):
	"""Pins the current output: --hic on simulated Hi-C pairs with chimeric read 1s."""
	g, idx, f1, f2, h1, h2 = pin_input
	res = H.run_map(idx, [h1, h2], ["--hic"])
	sam = res.sam()
	assert H.pair_invariant_errors(sam) == [] and not any(r.is_proper_pair for r in sam.records)
	H.check_pinned("c_hic.sam.gz", H.strip_pg(res.stdout))


def _orient_library(g, libs, n=50, n_unseed=6):
	"""Read 1 and read 2 records of libraries in the given orientations [(orient, insert mean)], n pairs
	each from anywhere in the genome plus n_unseed pairs whose read 1 or read 2 (alternately) has no
	seed, so only mate rescue in that orientation can place it."""
	r1, r2 = [], []
	for k, (orient, ins) in enumerate(libs):
		pairs = H.simulate_pairs(g, n + n_unseed, seed=130 + k, insert_mean=ins, insert_sd=30, sub_rate=0.01, indel_rate=0.001,
			qual=None, name_prefix=orient.lower())
		for i, p in enumerate(pairs):
			a, b, _ = C.reorient(p, orient)
			if i >= n and i & 1:
				a = C.substitute(a, C.unseedable_subs(len(a)))
			elif i >= n:
				b = C.substitute(b, C.unseedable_subs(len(b)))
			r1.append((p.name, a))
			r2.append((p.name, b))
	return r1, r2


ORIENT_RUNS = [
	("c_pe_orient.sam.gz", (("FR", 300), ("FF", 300), ("RR", 300), ("RF", 500)), 4),
	("c_pe_ff.sam.gz", (("FF", 300),), 1),
]


@pytest.mark.pinned
@pytest.mark.parametrize("name,libs,n_dir", ORIENT_RUNS, ids=["mixed", "ff"])
def test_orientation_output_pinned(pin_input, tmp_path, name, libs, n_dir):
	"""Pins the current output: paired-end SAM of a library mixing FR, FF, RR and RF pairs (statistics
	inferred for all four orientations), and of an FF-only library, both with mates that only rescue can
	place, so pairing, the rescue windows and the pair MAPQ are pinned in every orientation."""
	g, idx = pin_input[0], pin_input[1]
	r1, r2 = _orient_library(g, libs)
	f1 = H.write_fasta(str(tmp_path / "o1.fa"), r1)
	f2 = H.write_fasta(str(tmp_path / "o2.fa"), r2)
	res = H.run_map(idx, [f1, f2])
	cands, pct, bnd = C.stderr_pestat(res.stderr)
	assert len(cands) == 1 and len(bnd) == n_dir
	sam = res.sam()
	assert H.pair_invariant_errors(sam) == []
	assert sum(r.is_proper_pair for r in sam.records if r.is_primary) > 0
	H.check_pinned(name, H.strip_pg(res.stdout))
