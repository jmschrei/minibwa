"""Area B: presets, the short/long-read switch (-x, -l, --long, --adap), the
long-read heuristics (-W re-chaining, -g chain breaking, Z-drop splits,
inversions), --chain-only (hand-worked chains, and every simulated read's
best chain against a brute-force chaining DP over the anchors that
--dbg-anchor prints) and the alignment debug options. Reads are built from
the unique parts of helpers_b.BGenome so expected values can be worked out
by hand."""

import re

import pytest

import minibwa_helpers as H
import helpers_b as B
from helpers_b import bgen  # noqa: F401 (fixture)


def _map(bgen, tmp_path, reads, args=(), check=True):
	fa = H.write_fasta(str(tmp_path / "m.fa"), reads)
	return H.run_map(bgen.index, fa, list(args), check=check)


def _ins4_read(g, u, half):
	"""2*half bases from g[u:] with ACGT inserted in the middle (4 bp insertion)."""
	return g[u:u + half] + "ACGT" + g[u + half:u + 2 * half]


@pytest.mark.parametrize("args,long_mode", [
	((), True), (("-l", "405"), False), (("-l", "403"), True), (("--adap=no",), False),
	(("-x", "sr"), False), (("-x", "lr"), True), (("-x", "sr", "--long"), True), (("-x", "sr", "--long=yes"), True),
	(("-x", "sr", "--long", "--long=no"), False), (("--adap=no", "--long=y"), True), (("-x", "adap"), True),
])
def test_long_read_switch_404bp(bgen, tmp_path, args, long_mode):
	"""A 404 bp read with a 4 bp insertion. AS = 800 - min(12+8, 23+4) = 780 in every mode;
	ms charges the gap q+e = 14 in short-read mode (786) and q + e*log2(5) = 16.6 in
	long-read mode (783.4 -> 783). Reads longer than -l (325) are long in the adaptive
	mode (the default, -x adap); --adap=no, -x sr and a larger -l make them short; --long
	and -x lr make every read long."""
	g = bgen["b1"]
	rec = B.one_record(_map(bgen, tmp_path, [("L", _ins4_read(g, bgen.unique_start + 7000, 200))], args))
	assert rec.pos0 == bgen.unique_start + 7000 and rec.tags["AS"] == 780 and rec.tags["NM"] == 4, rec.line[:200]
	assert rec.tags["ms"] == (783 if long_mode else 786), (args, rec.line[:200])


@pytest.mark.parametrize("args,long_mode", [((), False), (("--long",), True), (("-x", "lr"), True), (("-l", "200"), True), (("-l", "205"), False)])
def test_long_read_switch_204bp(bgen, tmp_path, args, long_mode):
	"""A 204 bp read with a 4 bp insertion: AS 380; ms 386 short, 383 long (400 - 16.6)."""
	g = bgen["b1"]
	rec = B.one_record(_map(bgen, tmp_path, [("S", _ins4_read(g, bgen.unique_start + 7600, 100))], args))
	assert rec.tags["AS"] == 380 and rec.tags["ms"] == (383 if long_mode else 386), (args, rec.line[:200])


def test_preset_applies_before_other_options(bgen, tmp_path):
	"""-x is applied before every other option wherever it appears: '-w 5 -x sr' equals
	'-x sr -w 5', and both differ from '-x sr' (the narrow band splits a 40 bp deletion)."""
	g = bgen["b1"]
	u = bgen.unique_start + 8000
	while g[u + 74] == g[u + 114] or g[u + 75] == g[u + 115]:
		u += 1
	s = [("w", g[u:u + 75] + g[u + 115:u + 190])]
	a = H.strip_pg(_map(bgen, tmp_path, s, ("-w", "5", "-x", "sr")).stdout)
	b = H.strip_pg(_map(bgen, tmp_path, s, ("-x", "sr", "-w", "5")).stdout)
	c = H.strip_pg(_map(bgen, tmp_path, s, ("-x", "sr")).stdout)
	assert a == b != c
	assert len(H.sam_body(a).splitlines()) == 2 and len(H.sam_body(c).splitlines()) == 1


def test_preset_and_yes_no_errors(bgen, tmp_path):
	"""An unknown preset is an error (exit 1, message); --adap and --long take yes/no (or y/n)
	and warn on anything else without changing the setting."""
	r = _map(bgen, tmp_path, [("x", bgen["b1"][100:250])], ("-x", "nonsense"), check=False)
	assert r.returncode == 1 and "unknown preset 'nonsense'" in r.stderr and r.stdout == ""
	s = [("L", _ins4_read(bgen["b1"], bgen.unique_start + 7000, 200))]
	for opt in ("--adap=maybe", "--long=perhaps"):
		r = _map(bgen, tmp_path, s, (opt,))
		assert "only accepts 'yes' or 'no'" in r.stderr, opt
		assert B.one_record(r).tags["ms"] == 783, opt  # still the default adaptive mode


def test_chain_only_hand_worked(bgen, tmp_path):
	"""--chain-only writes PAF chains without a base alignment (no cg, NM, AS). Chain score
	s1 = first seed + each next seed's min(len, gap) - gap penalty; columns 10-11 from the seeds:
	a mismatch at 75: seeds 75 + 74 -> s1 149, cm 2, mlen 149, blen 150;
	a 30 bp deletion: 75 + 75 - (int)(0.152*30 + 0.5*log2 31) = 143, mlen 150, blen 180;
	a reverse-strand exact read: one 150 bp seed. MAPQ 40*(1-25/s1)*ln(s1) caps at 60."""
	g = bgen["b1"]
	u = bgen.unique_start + 9000
	m = g[u:u + 75] + B.flip_base(g[u + 75]) + g[u + 76:u + 150]
	v = u + 1000
	while g[v + 74] == g[v + 104] or g[v + 75] == g[v + 105]:
		v += 1
	d = g[v:v + 75] + g[v + 105:v + 180]
	w = u + 2000
	rc = H.revcomp(g[w:w + 150])
	L = len(g)
	want = {
		"m": "m\t150\t0\t150\t+\tb1\t%d\t%d\t%d\t149\t150\t60\ttp:A:P\ts1:i:149\tcm:i:2\ts2:i:0" % (L, u, u + 150),
		"d": "d\t150\t0\t150\t+\tb1\t%d\t%d\t%d\t150\t180\t60\ttp:A:P\ts1:i:143\tcm:i:2\ts2:i:0" % (L, v, v + 180),
		"rc": "rc\t150\t0\t150\t-\tb1\t%d\t%d\t%d\t150\t150\t60\ttp:A:P\ts1:i:150\tcm:i:1\ts2:i:0" % (L, w, w + 150),
	}
	reads = [("m", m), ("d", d), ("rc", rc)]
	out = _map(bgen, tmp_path, reads, ("--chain-only",)).stdout
	assert out.splitlines() == [want["m"], want["d"], want["rc"]]
	assert _map(bgen, tmp_path, reads, ("--chain-only", "-a")).stdout == out  # PAF is forced
	fa = H.write_fasta(str(tmp_path / "p1.fa"), [("m", m), ("d", d)])
	fb = H.write_fasta(str(tmp_path / "p2.fa"), [("m", H.revcomp(g[u + 300:u + 450])), ("d", H.revcomp(g[v + 300:v + 450]))])
	pe = H.run_map(bgen.index, [fa, fb], ["--chain-only"]).stdout.splitlines()  # no pairing either
	assert [l.split("\t")[0] for l in pe] == ["m/1", "m/2", "d/1", "d/2"]
	assert pe[0] == want["m"].replace("m\t", "m/1\t", 1) and pe[2] == want["d"].replace("d\t", "d/1\t", 1)


@pytest.fixture(scope="module")
def chain_runs(synth, tmp_path_factory):
	reads = H.simulate_reads(synth.genome, 250, length=150, seed=51, sub_rate=0.02, indel_rate=0.003, max_indel=6)
	fq = H.write_fastq(str(tmp_path_factory.mktemp("b_chain") / "c.fq"), reads)
	res = H.run_map(synth.index, fq, ["--chain-only", "--dbg-anchor"])
	return reads, res


def test_chain_only_matches_brute_force(synth, chain_runs):
	"""For every 150 bp read, the first (best) --chain-only line is the best chain of a
	brute-force chaining DP over the anchors --dbg-anchor prints (gap score 0.152*dd +
	0.5*log2(dd+1), band and distance from -w and -g adjusted for 150 bp: 107 and 150):
	same score; for reads from unique sequence (no equally good chain elsewhere) also the
	same anchor count, coordinates, strand and columns 10-11."""
	reads, res = chain_runs
	truth = dict((r.name, r) for r in reads)
	qlens = dict((r.name, len(r.seq)) for r in reads)
	ids = dict((n, i) for i, n in enumerate(synth.genome.names))
	anchors = B.parse_dbg_anchors(res.stderr, qlens, ids)
	first = {}
	for p in res.paf():
		first.setdefault(p.qname, p)
	bw, max_gap = B.adap_value(100, 500, 150), min(B.adap_value(100, 5000, 150), 150)
	assert (bw, max_gap) == (107, 150)
	n = n_full = 0
	for name, p in first.items():
		if p.is_unmapped:
			assert name not in anchors or B.chain_best(anchors[name], max_gap, bw)[0] < 25
			continue
		score, chain = B.chain_best(anchors[name], max_gap, bw)
		qs, qe, ts, te = B.chain_coords(chain, qlens[name])
		mlen, blen = B.chain_mlen_blen(chain)
		got = (p.tags["s1"], p.tags["cm"], p.qstart, p.qend, p.tstart, p.tend, p.strand, p.tname, p.nmatch, p.alen)
		want = (score, len(chain), qs, qe, ts, te, "-" if chain[0].rev else "+", chain[0].ctg, mlen, blen)
		t = truth[name]
		if synth.genome.is_unique(t.contig, t.pos, t.end, 200):
			assert got == want, (p.line, want)
			n_full += 1
		else:
			assert got[0] == want[0], (p.line, want)
		assert "cg" not in p.tags and "AS" not in p.tags
		n += 1
	assert n > 200 and n_full > 150


def test_dbg_anchors_are_exact_matches(synth, chain_runs):
	"""Every anchor --dbg-anchor prints is an exact match between the read (reverse-
	complemented for '-') and the genome."""
	reads, res = chain_runs
	seqs = dict((r.name, r.seq) for r in reads)
	n = 0
	for line in res.stderr.splitlines():
		f = line.split("\t")
		if f[0] != "AC":
			continue
		name, qs, strand, ctg, ts, length = f[1], int(f[2]), f[3], f[4], int(f[5]), int(f[6])
		q = seqs[name][qs:qs + length]
		if strand == "-":
			q = H.revcomp(q)
		assert synth.genome[ctg][ts:ts + length] == q, line
		n += 1
	assert n > 250


def test_g_breaks_chains(bgen, tmp_path):
	"""-x lr (-g 5000, -w 500, no length adjustment), --chain-only, a 5 kbp read whose
	middle 1000 bases are unrelated: two 2000 bp seeds 3000 apart on both sequences chain
	(2000 + 2000, no gap penalty; mlen 4000, blen 5000) unless -g 500 makes the distance too
	long. Re-chaining with the long band -W (20000) joins them again; -W 500 does not
	re-chain."""
	g = bgen["b2"]
	o = 500
	junk = H.Rng(91).dna(1000)
	junk = B.flip_base(g[o + 2000]) + junk[1:-1] + B.flip_base(g[o + 2999])
	s = [("j", g[o:o + 2000] + junk + g[o + 3000:o + 5000])]
	one = "j\t5000\t0\t5000\t+\tb2\t20000\t%d\t%d\t4000\t5000" % (o, o + 5000)
	for args in (("-W", "500"), ("-g", "500")):
		lines = _map(bgen, tmp_path, s, ("-x", "lr", "--chain-only") + args).stdout.splitlines()
		assert len(lines) == 1 and lines[0].startswith(one) and "\ts1:i:4000\tcm:i:2" in lines[0], (args, lines)
	lines = _map(bgen, tmp_path, s, ("-x", "lr", "--chain-only", "-W", "500", "-g", "500")).stdout.splitlines()
	got = sorted(tuple(l.split("\t")[2:4] + l.split("\t")[7:9]) for l in lines)
	assert got == sorted([("0", "2000", str(o), str(o + 2000)), ("3000", "5000", str(o + 3000), str(o + 5000))]), lines
	assert all("\ts1:i:2000\tcm:i:1\t" in l and "tp:A:P" in l for l in lines)


def test_rechaining_needs_W_above_twice_w(bgen, tmp_path):
	"""Re-chaining with the long band runs only when -W exceeds twice -w (500 with -x lr):
	with two 2000 bp seeds 2500 apart and -g 300 the initial chaining keeps them apart; -W 1001
	re-chains them into one chain (2000 + 2000; mlen 4000, blen 4500), -W 1000 does not. This
	threshold is minibwa's (mb_map_sai), not documented."""
	g = bgen["b2"]
	o = 500
	junk = H.Rng(91).dna(500)
	junk = B.flip_base(g[o + 2000]) + junk[1:-1] + B.flip_base(g[o + 2499])
	s = [("j", g[o:o + 2000] + junk + g[o + 2500:o + 4500])]
	lines = _map(bgen, tmp_path, s, ("-x", "lr", "--chain-only", "-g", "300", "-W", "1001")).stdout.splitlines()
	assert len(lines) == 1 and lines[0].startswith("j\t4500\t0\t4500\t+\tb2\t20000\t%d\t%d\t4000\t4500\t" % (o, o + 4500)) and "\ts1:i:4000\tcm:i:2" in lines[0], lines
	lines = _map(bgen, tmp_path, s, ("-x", "lr", "--chain-only", "-g", "300", "-W", "1000")).stdout.splitlines()
	assert sorted(tuple(l.split("\t")[2:4]) for l in lines) == [("0", "2000"), ("2500", "4500")], lines


def test_W_bridges_long_deletion(bgen, tmp_path):
	"""A 4 kbp read missing 2000 bp of the genome. Seeds on either side are too far off the
	diagonal for the initial band; re-chaining with -W re-joins them when -W exceeds the
	2000 bp shift: one alignment with a 2000 bp deletion (AS 8000 - (23 + 2000) = 5977),
	else a primary and a supplementary record."""
	g = bgen["b1"]
	o = 40000
	s = [("d", g[o:o + 2000] + g[o + 4000:o + 6000])]
	for args in ((), ("-W", "3000")):
		rec = B.one_record(_map(bgen, tmp_path, s, args))
		assert rec.pos0 <= o and re.match(r"^\d+M2000D\d+M$", rec.cigar) and rec.tags["AS"] == 5977 and rec.tags["NM"] == 2000, (args, rec.line[:200])
	recs = _map(bgen, tmp_path, s, ("-W", "1000")).sam().records
	assert sorted(r.flag & 0x900 for r in recs) == [0, 0x800]
	spans = sorted((r.pos0, r.ref_end0) for r in recs)
	assert abs(spans[0][0] - o) <= 2 and abs(spans[0][1] - (o + 2000)) <= 2 and abs(spans[1][0] - (o + 4000)) <= 2 and abs(spans[1][1] - (o + 6000)) <= 2, spans


def test_inversion(bgen, tmp_path):
	"""A 3.3 kbp read A + revcomp(B) + C from consecutive reference segments A (1500 bp),
	B (300 bp), C (1500 bp): A and C align forward, B in reverse (the inversion found
	between the two Z-drop-split halves), each as its own record with clips for the rest,
	all MAPQ 60 (an inversion takes the smaller MAPQ of its flanks). In PAF the inversion hit
	has no chain (s1 0, cm 0)."""
	g = bgen["b1"]
	o = 46100
	comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
	while comp[g[o + 1799]] == g[o + 1500] or comp[g[o + 1800]] == g[o + 1499]:  # no base fits both ways at a junction
		o += 1
	s = [("inv", g[o:o + 1500] + H.revcomp(g[o + 1500:o + 1800]) + g[o + 1800:o + 3300])]
	recs = _map(bgen, tmp_path, s).sam().records
	got = sorted((r.strand, r.pos0, r.ref_end0, r.cigar.replace("H", "S"), r.mapq) for r in recs)
	assert got == sorted([("+", o, o + 1500, "1500M1800S", 60), ("+", o + 1800, o + 3300, "1800S1500M", 60),
		("-", o + 1500, o + 1800, "1500S300M1500S", 60)]), [r.line[:120] for r in recs]
	assert sum(r.is_primary for r in recs) == 1 and all("SA" in r.tags for r in recs)
	paf = _map(bgen, tmp_path, s, ("-f",)).paf()
	inv = [p for p in paf if p.strand == "-"]
	assert len(paf) == 3 and len(inv) == 1 and (inv[0].qstart, inv[0].qend, inv[0].tags["s1"], inv[0].tags["cm"]) == (1500, 1800, 0, 0)


def test_zdrop_splits_unrelated_middle(bgen, tmp_path):
	"""A 3 kbp read with its middle 1000 bases replaced by unrelated sequence: the two flanks
	chain, the alignment through the middle drops (Z-drop) and the chain is split, giving
	one primary and one supplementary record of 1000 matches each."""
	g = bgen["b2"]
	o = 5600
	junk = H.Rng(92).dna(1000)
	junk = B.flip_base(g[o + 1000]) + junk[1:-1] + B.flip_base(g[o + 1999])
	s = [("z", g[o:o + 1000] + junk + g[o + 2000:o + 3000])]
	res = _map(bgen, tmp_path, s, ("--dbg-an-pos",))
	recs = res.sam().records
	got = sorted((r.pos0, r.cigar.replace("H", "S")) for r in recs)
	assert got == [(o, "1000M2000S"), (o + 2000, "2000S1000M")] and sorted(r.flag & 0x900 for r in recs) == [0, 0x800], got
	ad = [l.split("\t") for l in res.stderr.splitlines() if l.startswith("AD\t")]
	assert any(f[7] == "1" for f in ad), ad  # the gap fill was Z-dropped


def test_dbg_an_pos_and_aln_seq(bgen, tmp_path):
	"""--dbg-an-pos prints the anchors of each hit (AF: hit, contig, last target base, last
	query base, gap to the previous anchor, length; all exact matches) and the gap fills (AD);
	--dbg-aln-seq prints each DP's parameters, target and query (reversed for a left
	extension) and result."""
	g = bgen["b1"]
	u = bgen.unique_start + 12000
	while g[u + 99] == g[u + 102] or g[u + 100] == g[u + 103]:
		u += 1
	# a mismatch at 40 (two seeds in one ungapped stretch) and a 3 bp deletion after base 100
	read = g[u:u + 40] + B.flip_base(g[u + 40]) + g[u + 41:u + 100] + g[u + 103:u + 153]
	res = _map(bgen, tmp_path, [("x", read)], ("--dbg-an-pos", "--dbg-aln-seq", "-O", "6", "-E", "3"))
	af = [l.split("\t") for l in res.stderr.splitlines() if l.startswith("AF\t")]
	assert len(af) >= 2
	prev = None
	for f in af:
		ctg, tpos, qpos, gap, length = f[2], int(f[3]), int(f[4]), int(f[5]), int(f[6])
		assert ctg == "b1" and g[tpos + 1 - length:tpos + 1] == read[qpos + 1 - length:qpos + 1]
		assert gap == (0 if prev is None else (qpos - prev[1]) - (tpos - prev[0]))
		prev = (tpos, qpos)
	lines = res.stderr.splitlines()
	heads = [i for i, l in enumerate(lines) if l.startswith("===> ")]
	assert heads
	for i in heads:
		assert lines[i].startswith("===> q=(6,6), e=(3,3), bw=")
		t, q, sc = lines[i + 1], lines[i + 2], lines[i + 3]
		assert t in g or t[::-1] in g
		assert q in read or q[::-1] in read
		assert re.match(r"^score=-?\d+, max=\d+, cigar=(\d+[MIDN])*$", sc), sc
	rec = B.one_record(res)
	assert (rec.pos0, rec.cigar, rec.tags["AS"]) == (u, "100M3D50M", 298 - 8 - 15)
