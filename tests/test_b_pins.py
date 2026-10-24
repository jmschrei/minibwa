"""Area B byte pins: `minibwa map` output on fixed simulated inputs, for the
code paths a speed optimization of chaining and alignment is most likely to
touch (chaining, extension and gap filling, Z-drop and splits, secondary
selection, MAPQ, every base-alignment tag, the presets and the scoring
options). These pin the current output; the other test_b_* files check the
same paths against independently derived values."""

import pytest

import minibwa_helpers as H


@pytest.fixture(scope="module")
def pin_inputs(synth, tmp_path_factory):
	"""Short reads (50-320 bp, substitutions, indels up to 10 bp, N), chimeric reads (two
	loci), accurate long reads (330 bp - 6 kbp) and noisier long reads (0.8-3 kbp, 3%
	substitutions and 2% indels up to 30 bp, which exercise the long-read end clean-up in
	mm_fix_bad_ends) on the standard synthetic genome."""
	g = synth.genome
	d = tmp_path_factory.mktemp("b_pins")
	se = H.simulate_reads(g, 300, length=(50, 320), seed=61, sub_rate=0.02, indel_rate=0.004, max_indel=10, n_rate=0.002)
	chim = [p.r1 for p in H.simulate_hic_pairs(g, 40, read_length=150, seed=63, chimeric_frac=1.0, unique=False, name_prefix="c")]
	lr = H.simulate_long_reads(g, 40, length=(330, 6000), seed=62, sub_rate=0.01, indel_rate=0.004, max_indel=12)
	noisy = H.simulate_long_reads(g, 200, length=(800, 3000), seed=77, sub_rate=0.03, indel_rate=0.02, max_indel=30)
	return dict(se=H.write_fastq(str(d / "se.fq"), se + chim), lr=H.write_fastq(str(d / "lr.fq"), lr + _structural_reads(g)),
		ends=H.write_fastq(str(d / "ends.fq"), _end_indel_reads(g)), noisy=H.write_fastq(str(d / "noisy.fq"), noisy))


def _end_indel_reads(g):
	"""150 bp reads with a 5-60 bp deletion or a 1-12 bp insertion 8-45 bases from one end:
	where the extension window, the end bonus and the gap costs decide between clipping and
	aligning through the indel."""
	rng = H.Rng(65)
	out = []
	for i in range(120):
		c = g.names[i % 2]
		seq = g[c]
		o = rng.below(len(seq) - 400)
		tail = 8 + rng.below(38)
		head = 150 - tail if i % 4 < 2 else tail  # the indel near the right or the left end
		if i % 3 == 2:  # insertion
			m = 1 + rng.below(12)
			r = seq[o:o + head - m] + rng.dna(m) + seq[o + head - m:o + 150 - m]
		else:  # deletion
			m = 5 + rng.below(56)
			r = seq[o:o + head] + seq[o + head + m:o + 150 + m]
		if i % 2:
			r = H.revcomp(r)  # mapped to the reverse strand
		out.append(("e%d" % i, r, "I" * len(r)))
	return out


def _structural_reads(g):
	"""Long reads with an inversion, a long deletion, an unrelated middle or two loci, so the
	pins cover re-chaining, Z-drop splits, inversion alignment and supplementary records."""
	rng = H.Rng(64)
	out = []
	for i in range(24):
		kind = i % 4
		c = g.names[i % 2]
		seq = g[c]
		o = rng.below(len(seq) - 9000)
		a, b = 1000 + rng.below(1000), 1000 + rng.below(1000)
		if kind == 0:  # inversion of 250-600 bp
			m = 250 + rng.below(350)
			r = seq[o:o + a] + H.revcomp(seq[o + a:o + a + m]) + seq[o + a + m:o + a + m + b]
		elif kind == 1:  # deletion of 500-3000 bp
			m = 500 + rng.below(2500)
			r = seq[o:o + a] + seq[o + a + m:o + a + m + b]
		elif kind == 2:  # 300-1000 unrelated bases in place of the reference
			m = 300 + rng.below(700)
			r = seq[o:o + a] + rng.dna(m) + seq[o + a + m:o + a + m + b]
		else:  # two loci on different contigs, the second reverse-complemented
			c2 = g.names[(i + 1) % 2]
			o2 = rng.below(len(g[c2]) - 3000)
			r = seq[o:o + a] + H.revcomp(g[c2][o2:o2 + b])
		r = H.mutate(rng, r, 0.003)
		out.append(("sv%d" % i, r, "I" * len(r)))
	return out


PINS = [
	("b_se_cs.sam.gz", "se", ["-b", "cs", "--outn", "5"]),
	("b_se_eqx_ds.paf.gz", "se", ["-f", "--eqx", "-b", "ds", "--outn", "5"]),
	("b_se_scoring.sam.gz", "se", ["-A", "3", "-B", "5", "-O", "6,30", "-E", "3,1", "-w", "30", "-m", "30", "-s", "40", "-p", "0.8", "-N", "3", "--outn", "3", "-b", "MD"]),
	("b_se_long_mode.paf.gz", "se", ["--long", "-f", "-b", "cs"]),
	("b_lr_md.sam.gz", "lr", ["-b", "MD", "--outn", "2"]),
	("b_lr_preset_lr.paf.gz", "lr", ["-x", "lr", "-f", "-b", "cs"]),
	("b_lr_preset_sr.sam.gz", "lr", ["-x", "sr", "-b", "ds"]),
	("b_lr_bands.sam.gz", "lr", ["-W", "2000", "-g", "2000", "-l", "1000", "-w", "50", "--eqx"]),
	("b_se_chain_only.paf.gz", "se", ["--chain-only", "--outn", "5"]),
	("b_end_indels.sam.gz", "ends", ["-b", "cs"]),
	("b_end_indels_scores.paf.gz", "ends", ["-f", "-O", "6,20", "-E", "3,1", "-B", "6", "-w", "20"]),
	("b_lr_chain_only.paf.gz", "lr", ["--chain-only"]),
	("b_lr_noisy.paf.gz", "noisy", ["-f"]),
]


@pytest.mark.pinned
@pytest.mark.parametrize("pin,inp,args", PINS, ids=[p[0] for p in PINS])
def test_pinned_output(synth, pin_inputs, pin, inp, args):
	"""Pins the current output of `minibwa map` (without the @PG line) on fixed inputs."""
	res = H.run_map(synth.index, pin_inputs[inp], args)
	H.check_pinned(pin, H.strip_pg(res.stdout))
