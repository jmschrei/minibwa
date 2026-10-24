"""Area A: SMEMs, seed intervals and anchors against brute force.

Expected values come from helpers_a: SMEMs found by extending every read
position as far as its occurrence count allows (counts from a naive suffix
array of the indexed text), sub-SMEMs by the rule minibwa documents in
seed.c, and anchors from the SA values of each interval mapped to contigs.
The checks run through the C driver (mb_bwt_smem, mb_bwt_smem_batch,
mb_seed_intv, mb_seed_intv_batch, mb_anchor), through `minibwa fastmap`, and
through `minibwa map --dbg-seed --dbg-anchor`, which prints the intervals and
anchors of the mapping path itself.
"""

import os

import pytest

import minibwa_helpers as H
import helpers_a as A


@pytest.fixture(scope="module")
def tiny(index_of, tmp_path_factory):
	return A.Tiny(index_of, str(tmp_path_factory.mktemp("tiny")))


# ---------------------------------------------------------------------------
# mb_bwt_smem and mb_bwt_smem_batch
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("min_len,min_occ,cache", [(19, 1, 10), (10, 1, 0), (5, 1, 10), (12, 3, 10), (25, 2, 10)])
def test_smem_equals_brute_force(tiny, min_len, min_occ, cache):
	"""mb_bwt_smem finds exactly the SMEMs of each read (maximal matches occurring at least
	min_occ times, at least min_len long, not contained in another), with the SA interval of
	the match and of its reverse complement, over the whole read and over its middle half;
	N and IUPAC bases in reads break matches. With and without the 10-mer cache. The driver
	checks that mb_bwt_smem_batch finds the same sets in batches of 1, 3, 50 and 1000
	entries, with whole reads and middle halves mixed in one batch."""
	r = H.assert_c_case("drv_a_bwt", "smem", tiny.index + ".mbw", tiny.reads_fa, min_len, min_occ, cache, 1, 3, 50, 1000)
	whole = A.intervals_by_read(r.stdout, "M")
	mid = A.intervals_by_read(r.stdout, "m")
	n = 0
	for name, _ in tiny.reads:
		q = tiny.codes[name]
		want = set(A.smems(tiny.ref, q, min_len, min_occ))
		assert whole.get(name, set()) == want, name
		st, en = len(q) // 4, len(q) - len(q) // 4
		assert mid.get(name, set()) == set(A.smems(tiny.ref, q, min_len, min_occ, st, en)), name
		n += len(want)
	assert n > 50


# ---------------------------------------------------------------------------
# mb_seed_intv, mb_seed_intv_batch, mb_anchor
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("min_len,max_sub_occ", [(19, 10), (12, 10), (15, 20), (19, 0)])
def test_seed_intervals_equal_brute_force(tiny, min_len, max_sub_occ):
	"""mb_seed_intv returns the SMEMs plus, inside each SMEM at least 2*min_len long that
	occurs at most max_sub_occ times, the SMEMs that occur more often than it and are at
	least max(len/2, min_len) long. mb_seed_intv_batch returns the same set for reads
	processed in groups of 1, 7, 60 and all (checked in C)."""
	r = H.assert_c_case("drv_a_bwt", "seed", tiny.index, tiny.reads_fa, min_len, max_sub_occ, 1, 7, 60, 1000)
	got = A.intervals_by_read(r.stdout, "D")
	n_sub = 0
	for name, _ in tiny.reads:
		q = tiny.codes[name]
		want = A.seed_intervals(tiny.ref, q, min_len, max_sub_occ)
		assert got.get(name, set()) == want, name
		n_sub += len(want - set(A.smems(tiny.ref, q, min_len, 1)))
	if max_sub_occ > 0:
		assert n_sub > 10  # sub-SMEMs were exercised
	else:
		assert n_sub == 0


def _anchor_rows(text):
	got = {}
	for f in A.parse_rows(text, "A"):
		got.setdefault(f[0], set()).add((int(f[1]), "+-".index(f[2]), int(f[3]), int(f[4]), int(f[5])))
	return got


@pytest.mark.parametrize("min_len,max_occ", [(19, 250), (12, 5), (10, 3), (15, 1)])
def test_anchors_equal_brute_force(tiny, min_len, max_occ):
	"""mb_anchor turns every interval into anchors at its SA positions (all of them, or
	max_occ evenly spaced rows when it occurs more often), drops occurrences spanning two
	contigs or the strand boundary, and drops an anchor contained in another on the same
	contig, strand and diagonal; it returns 1.0 and sorts by contig, strand and end. The
	driver checks the order and that mb_anchor_sort restores it from a shuffled copy."""
	r = H.assert_c_case("drv_a_bwt", "anchor", tiny.index, tiny.reads_fa, min_len, 10, max_occ)
	got = _anchor_rows(r.stdout)
	n_sampled = 0
	for name, _ in tiny.reads:
		q = tiny.codes[name]
		seeds = A.seed_intervals(tiny.ref, q, min_len, 10)
		assert got.get(name, set()) == A.anchors(tiny.ref, len(q), seeds, max_occ), name
		n_sampled += sum(1 for s in seeds if s[4] > max_occ)
	if max_occ < 10:
		assert n_sampled > 0


def test_anchor_corner_cases(tiny):
	"""Every anchor lies inside one contig and matches the text there exactly, also for
	reads crossing a contig junction or the strand boundary of the text; reads through the reference N run match the random
	bases the index stores there; reads with N, IUPAC codes or '.' give no anchor across
	those bases; an exact copy of a contig gives one anchor of its full length."""
	r = H.assert_c_case("drv_a_bwt", "anchor", tiny.index, tiny.reads_fa, 19, 10, 250)
	got = _anchor_rows(r.stdout)
	g = dict(tiny.contigs)
	seqs = dict(tiny.reads)
	for name, seq in tiny.reads:
		for (tid, rev, qs, ts, ln) in got.get(name, ()):
			ctg = tiny.ref.names[tid]
			assert 0 <= ts and ts + ln <= len(g[ctg])
			part = seq[qs:qs + ln].upper()
			assert "N" not in part and set(part) <= set("ACGT")
			text = tiny.ref.fwd[tiny.ref.offs[tid] + ts:tiny.ref.offs[tid] + ts + ln]
			assert (part if not rev else H.revcomp(part)) == text, name
	s2 = [n for n in seqs if n.startswith("exact_s2")][0]
	assert (3, 1, 0, 0, 33) in got[s2]
	nrun = [n for n in seqs if n.startswith("nrun_")]
	assert all(n in got for n in nrun)
	# the reads taken from the text cross t1's 30 bp N run in one anchor
	assert any(ts <= 1500 and ts + ln >= 1530 for n in nrun[:2] for (_, _, _, ts, ln) in got[n])


# ---------------------------------------------------------------------------
# minibwa fastmap
# ---------------------------------------------------------------------------


def _fastmap_expect(ref, codes, l, s, w):
	out = []
	for (qs, qe, x0, x1, size) in A.smems(ref, codes, l, s):
		if size > w:
			pos = None
		else:
			pos = []
			for j in range(size):
				p = ref.sa_row(x0 + j)
				c = ref.intv2cid(p, p + qe - qs)
				pos.append("%s:%s%d" % (ref.names[c[0]], "+-"[c[2]], c[1] + 1) if isinstance(c, tuple) else ".")
		out.append((qs, qe, size, pos))
	return out


@pytest.mark.parametrize("args,l,s,w", [
	([], 19, 1, 20),
	(["-l", "10", "-s", "2", "-w", "5", "-b", "7"], 10, 2, 5),
	(["-l", "5", "-w", "1000", "-b", "1000"], 5, 1, 1000),
	(["-l", "12", "-s", "3", "-w", "0"], 12, 3, 0),
])
def test_fastmap_equals_brute_force(tiny, args, l, s, w):
	"""minibwa fastmap prints, for each read, every SMEM (start, end, occurrence count) and,
	when it occurs at most -w times, each occurrence in suffix-array order as
	contig:strand+1-based start, or '.' when it spans two contigs or the strand boundary;
	-b > 1 (mb_bwt_smem_batch) prints the same as -b 1 (mb_bwt_smem)."""
	r = H.run_minibwa(["fastmap"] + args + [tiny.index, tiny.reads_fa], check=True)
	got = A.parse_fastmap(r.stdout)
	assert [x[0] for x in got] == [n for n, _ in tiny.reads]
	n_dot = 0
	for (name, ln, ems), (_, seq) in zip(got, tiny.reads):
		assert ln == len(seq)
		want = _fastmap_expect(tiny.ref, tiny.codes[name], l, s, w)
		assert ems == want, name
		n_dot += sum(p.count(".") for _, _, _, p in ems if p)
	if w >= 20 and l <= 19:
		assert n_dot > 0  # a junction read was exercised


def test_fastmap_stdin_gzip_and_usage(tiny, tmp_path):
	"""fastmap reads gzip and stdin ('-') like plain files; --help prints its usage."""
	plain = H.run_minibwa(["fastmap", "-l", "15", tiny.index, tiny.reads_fa], check=True).stdout
	gz = H.write_fasta(str(tmp_path / "r.fa.gz"), tiny.reads)
	assert H.run_minibwa(["fastmap", "-l", "15", tiny.index, gz], check=True).stdout == plain
	with open(tiny.reads_fa, "rb") as f:
		data = f.read()
	assert H.run_minibwa(["fastmap", "-l", "15", "-b", "4", tiny.index, "-"], stdin=data, check=True).stdout == plain
	r = H.run_minibwa(["fastmap", "--help"])
	assert r.returncode == 0 and "Usage: minibwa fastmap" in r.stdout and "-w INT" in r.stdout
	r = H.run_minibwa(["fastmap", tiny.index])
	assert r.returncode == 1 and "Usage: minibwa fastmap" in r.stderr


# ---------------------------------------------------------------------------
# The mapping path: --dbg-seed and --dbg-anchor
# ---------------------------------------------------------------------------


def _dbg(stderr, ref):
	sd, ac = {}, {}
	for f in A.dbg_lines(stderr, "SD"):
		sd.setdefault(f[0], set()).add((int(f[1]), int(f[2]), int(f[3])))
	for f in A.dbg_lines(stderr, "AC"):
		ac.setdefault(f[0], set()).add((ref.names.index(f[3]), "+-".index(f[2]), int(f[1]), int(f[4]), int(f[5])))
	return sd, ac


@pytest.mark.parametrize("args,k,c", [
	([], 19, 250),
	(["-k", "12", "-c", "3"], 12, 3),
	(["-k25", "-c", "0.005k"], 25, 5),
	(["-x", "sr", "-k", "15", "-c", "1"], 15, 1),
])
def test_map_seeds_and_anchors(tiny, threads, args, k, c):
	"""minibwa map --dbg-seed prints, per read, the seed intervals of the mapping path
	(start, end, occurrence count): the SMEMs of length >= -k plus the sub-SMEMs (at most 10
	occurrences); --dbg-anchor prints the anchors after -c sampling and deduplication. Both
	equal the brute-force sets, for any -t; -c takes k/m suffixes."""
	r = H.run_map(tiny.index, tiny.reads_fa, args + ["--dbg-seed", "--dbg-anchor"])
	sd, ac = _dbg(r.stderr, tiny.ref)
	for name, _ in tiny.reads:
		q = tiny.codes[name]
		seeds = A.seed_intervals(tiny.ref, q, k, 10)
		assert sd.get(name, set()) == set((s, e, n) for (s, e, _, _, n) in seeds), name
		assert ac.get(name, set()) == A.anchors(tiny.ref, len(q), seeds, c), name
	assert len(r.sam().records) >= len(tiny.reads)


def test_map_read_encoding(tiny, tmp_path):
	"""The mapping path encodes reads with A/C/G/T in either case as bases and anything else
	as N: a lowercase read seeds like its uppercase form, U/u seed like T/t (bseq converts
	them), and an IUPAC code or '.' breaks seeds exactly as N does. SEQ keeps the case."""
	base = [(n, s) for n, s in tiny.reads if n.startswith("sim")][:12]
	variants = []
	for n, s in base:
		variants.append((n + "_up", s.upper()))
		variants.append((n + "_lo", s.lower()))
		variants.append((n + "_u", s.upper().replace("T", "U")))
		amb = list(s.upper())
		for p in (len(s) // 3, len(s) // 2):
			amb[p] = "R"
		variants.append((n + "_iupac", "".join(amb)))
		amb[len(s) // 3] = "."
		amb[len(s) // 2] = "N"
		variants.append((n + "_n", "".join(amb)))
	fa = A.write_reads_fasta(str(tmp_path / "v.fa"), variants)
	r = H.run_map(tiny.index, fa, ["--dbg-seed"])
	sd, _ = _dbg(r.stderr, tiny.ref)
	n_seeded = 0
	for n, _ in base:
		up = sd.get(n + "_up", set())
		assert sd.get(n + "_lo", set()) == up and sd.get(n + "_u", set()) == up, n
		n_seeded += len(up) > 0
		assert sd.get(n + "_iupac") == sd.get(n + "_n"), n
		want = A.seed_intervals(tiny.ref, A.nt4_codes(dict(variants)[n + "_n"]), 19, 10)
		assert sd.get(n + "_n", set()) == set((s, e, c) for (s, e, _, _, c) in want), n
	assert n_seeded >= 8
	recs = dict((x.qname, x) for x in r.sam().records if x.is_primary)
	for n, s in variants:
		rec = recs[n]
		want = s.replace("U", "T").replace("u", "t")
		assert rec.seq == (H.revcomp(want) if rec.is_reverse else want), n


def test_map_seeds_real_reads(chrm, tmp_path):
	"""On 150 real chrM read pairs (mapped as 300 single reads), the mapping path's seed
	intervals and anchors equal the brute-force sets computed on the chrM text."""
	ref = A.RefText(list(chrm.genome.items()))
	recs = []
	for fn, tag in ((chrm.read1, "/1"), (chrm.read2, "/2")):
		recs += [(n + tag, s) for n, s, _, _ in H.read_fastx(fn)[:150]]
	fa = A.write_reads_fasta(str(tmp_path / "r.fa"), recs)
	r = H.run_map(chrm.index, fa, ["--pe=no", "--dbg-seed", "--dbg-anchor"])
	sd, ac = _dbg(r.stderr, ref)
	for name, seq in recs:
		q = A.nt4_codes(seq)
		seeds = A.seed_intervals(ref, q, 19, 10)
		assert sd.get(name, set()) == set((s, e, c) for (s, e, _, _, c) in seeds), name
		assert ac.get(name, set()) == A.anchors(ref, len(q), seeds, 250), name
	assert len(sd) >= 290
