"""Area A: the index lookups seeding is built on, against a naive suffix array.

The tiny genome (helpers_a.make_tiny) is small enough to sort every suffix of
the indexed text T = F + revcomp(F) in Python, so ranks, SA values, SA
intervals, the 10-mer cache and the .l2b contig lookups are checked against
values derived from the FASTA alone. Covers bwt.c lookups (rank, extension,
cache, SA), l2bit.c lookups, kom_nt4_table, and the `getref` and `bench`
subcommands.
"""

import os
import re
import shutil

import pytest

import minibwa_helpers as H
import helpers_a as A


@pytest.fixture(scope="module")
def tiny(index_of, tmp_path_factory):
	return A.Tiny(index_of, str(tmp_path_factory.mktemp("tiny")))


@pytest.fixture(scope="module")
def edge(index_of, tmp_path_factory):
	contigs = A.make_edge()
	fa = A.write_contigs(str(tmp_path_factory.mktemp("edge") / "edge.fa"), contigs, line_width=7)
	return contigs, index_of(fa), A.RefText(contigs, with_sa=False)


# ---------------------------------------------------------------------------
# BWT: occurrence counts, SA, bidirectional extension, k-mer cache
# ---------------------------------------------------------------------------


def test_rank_equals_naive_occurrence_counts(tiny):
	"""mb_bwt_rank1a(k) equals the A/C/G/T counts in BWT rows [0, k) of the naive BWT for
	every k in [0, n+1]; the header (length, primary row, cumulative counts) matches. The
	driver also checks mb_bwt_rank11 at every k and mb_bwt_rank2a for all pairs within 300
	rows and 20,000 random pairs (k <= l)."""
	ref = tiny.ref
	r = H.assert_c_case("drv_a_bwt", "rank", tiny.index + ".mbw")
	head = A.parse_rows(r.stdout, "N")
	assert head == [[str(ref.n), str(ref.primary())] + [str(x) for x in ref.l2()[1:]]]
	occ = ref.occ_table()
	rows = A.parse_rows(r.stdout, "R")
	assert len(rows) == ref.n + 2
	for f in rows:
		k = int(f[0])
		assert tuple(int(x) for x in f[1:]) == occ[k], "rank at row %d" % k


def test_sa_equals_naive_suffix_array(tiny):
	"""mb_bwt_sa(row) is the start of the row's suffix in T for every row (row 0, the
	sentinel, may read as -1 or n); the driver checks mb_bwt_sa_batch against it in
	chunks of 1, 3, 64 and 1000 rows and on shuffled rows with repeats."""
	ref = tiny.ref
	r = H.assert_c_case("drv_a_bwt", "sa", tiny.index + ".mbw", 1, 3, 64, 1000)
	rows = A.parse_rows(r.stdout, "S")
	assert len(rows) == ref.n + 1
	assert int(rows[0][1]) in (-1, ref.n)
	got = [int(f[1]) for f in rows[1:]]
	assert got == ref.sa


def test_extension_of_every_short_pattern(tiny):
	"""For all 1,364 patterns of 1-5 bases, the interval built by backward extension
	(mb_bwt_set_intv then mb_bwt_extend with is_back=1) has the size of the pattern's
	occurrence count in T, starts at its naive SA interval and carries the interval of
	its reverse complement; forward extension gives the same interval (checked in C)."""
	ref = tiny.ref
	r = H.assert_c_case("drv_a_bwt", "extend", tiny.index + ".mbw", 5)
	rows = A.parse_rows(r.stdout, "X")
	assert len(rows) == 4 + 16 + 64 + 256 + 1024
	n_absent = 0
	for p, x0, x1, size in rows:
		w = ref.bi_interval(p)
		assert int(size) == w[2], p
		if w[2] > 0:
			assert (int(x0), int(x1)) == w[:2], p
		else:
			n_absent += 1
	assert 0 < n_absent < 1024  # the test needs both kinds


def test_kmer_cache(tiny):
	"""mb_bwt_cache(4) stores, at index sum(base_i << 2i), the interval of every 4-mer as in
	the naive suffix array. mb_bwt_cache(10), the cache minibwa map builds at index load,
	agrees with base-by-base backward search for all 4^10 10-mers (checked in C), and its
	occurrence counts sum to the number of 10-mer positions in T."""
	ref = tiny.ref
	r = H.assert_c_case("drv_a_bwt", "cache", tiny.index + ".mbw", 4, 1)
	rows = A.parse_rows(r.stdout, "K")
	assert len(rows) == 256
	for z, x0, x1, size in rows:
		z = int(z)
		kmer = "".join("ACGT"[z >> 2 * i & 3] for i in range(4))
		w = ref.bi_interval(kmer)
		assert int(size) == w[2], kmer
		if w[2]:
			assert (int(x0), int(x1)) == w[:2], kmer
	assert A.parse_rows(r.stdout, "O") == [[str(ref.n - 3)]]
	r = H.assert_c_case("drv_a_bwt", "cache", tiny.index + ".mbw", 10, 0)
	assert A.parse_rows(r.stdout, "O") == [[str(ref.n - 9)]]


def test_index_loading_and_failures(tiny, tmp_path):
	"""mb_bwt_load and l2b_load read a good index and return NULL for a missing file, a file
	of the other kind or a file of junk, and so does mb_idx_load for a missing prefix;
	minibwa map aborts with an error when the .mbw is missing or is not an .mbw."""
	bad = str(tmp_path / "bad.bin")
	with open(bad, "wb") as f:
		f.write(b"not an index at all\n")
	H.assert_c_case("drv_a_bwt", "load", tiny.index + ".mbw", tiny.index + ".l2b", bad)
	prefix = str(tmp_path / "copy")
	shutil.copy(tiny.index + ".l2b", prefix + ".l2b")
	reads = A.write_reads_fasta(str(tmp_path / "r.fa"), tiny.reads[:3])
	for mbw in (None, tiny.index + ".l2b"):  # missing, then the wrong magic number
		if mbw:
			shutil.copy(mbw, prefix + ".mbw")
		r = A.run_minibwa_abort(["map", prefix, reads])
		assert r.returncode != 0 and "failed to load the index" in r.stderr
		assert r.stdout == ""


# ---------------------------------------------------------------------------
# .l2b lookups
# ---------------------------------------------------------------------------


def _check_l2b(out, contigs, ref, comments):
	head = A.parse_rows(out, "H")
	assert head == [[str(len(contigs)), str(ref.tot), str(len(ref.ambi))]]
	ctg = A.parse_rows(out, "C")
	assert len(ctg) == len(contigs)
	for i, f in enumerate(ctg):
		name, seq = contigs[i]
		comm = "\t".join(f[4:])
		assert f[:4] == [str(i), name, str(len(seq)), str(ref.offs[i])]
		assert comm == comments.get(name, "*")
	assert [tuple(int(x) for x in f) for f in A.parse_rows(out, "B")] == ref.ambi
	seqs = dict((int(f[0]), f[1] if len(f) > 1 else "") for f in A.parse_rows(out, "G"))
	for i in range(len(contigs)):
		assert seqs[i] == "".join(str(c) for c in ref.contig_codes(i)), contigs[i][0]
	# l2b_getambi returns the index range of the ambiguous intervals overlapping a window
	n_win = 0
	for f in A.parse_rows(out, "W"):
		tid, st, en, ret, n = (int(x) for x in f)
		lo, hi = ref.offs[tid] + st, ref.offs[tid] + min(en, ref.lens[tid])
		want = [k for k, (a, b) in enumerate(ref.ambi) if a < hi and lo < b] if lo < hi else []
		if want:
			assert (ret, n) == (want[0], len(want)), f
			assert want == list(range(want[0], want[0] + len(want)))
		else:
			assert (ret, n) == (-1, 0), f
		n_win += 1
	assert n_win > 0
	# l2b_intv2cid on every start of both strands
	n_ok = 0
	for f in A.parse_rows(out, "I"):
		st, en, ret = int(f[0]), int(f[1]), int(f[2])
		want = ref.intv2cid(st, en)
		if isinstance(want, tuple):
			assert (ret, int(f[3]), int(f[4])) == want, f
			n_ok += 1
		else:
			assert ret == want, f
	return n_ok


def test_l2b_lookups_tiny(tiny):
	"""l2b_load: contig names, comments, lengths, offsets and the ambiguous intervals
	(maximal runs of non-ACGT bases, IUPAC codes included); l2b_getseq gives 0-3 for
	A/C/G/T in either case and 4 for ambiguous bases, clamps to the contig and agrees on
	sub-windows; l2b_getambi finds the overlapping intervals; l2b_intv2cid maps every
	text interval of 1, 3 and 40 bases to its contig, 0-based start and strand, or to -1
	(spans two contigs), -2 (spans the strand boundary) or -3 (past the text)."""
	r = H.assert_c_case("drv_a_bwt", "l2b", tiny.index + ".l2b", 1, 3, 40)
	n_ok = _check_l2b(r.stdout, tiny.contigs, tiny.ref, A.TINY_COMMENTS)
	assert n_ok > 3 * tiny.ref.tot


def test_l2b_lookups_edge(edge):
	"""The same lookups on contigs of 0, 1, 31, 32, 33, 64 and 100 bp."""
	contigs, prefix, ref = edge
	r = H.assert_c_case("drv_a_bwt", "l2b", prefix + ".l2b", 1, 2, 31, 33)
	assert _check_l2b(r.stdout, contigs, ref, {}) > 0


def test_nt4_table():
	"""kom_nt4_table: A/C/G/T in either case -> 0-3, the codes 0-3 themselves -> 0-3 (the
	API accepts 0/1/2/3-encoded queries), every other byte -> 4."""
	r = H.assert_c_case("drv_a_bwt", "nt4")
	table = dict((int(a), int(b)) for a, b in A.parse_rows(r.stdout, "T"))
	assert len(table) == 256
	for i in range(256):
		c = chr(i)
		want = "ACGT".index(c.upper()) if c in "ACGTacgt" else i if i < 4 else 4
		assert table[i] == want, i


# ---------------------------------------------------------------------------
# getref
# ---------------------------------------------------------------------------


def test_getref_tiny(tiny):
	"""minibwa getref prints every contig on one line, uppercase, with N for every
	ambiguous base (IUPAC codes and lowercase n included)."""
	r = H.run_minibwa(["getref", tiny.index + ".l2b"], check=True)
	want = "".join(">%s\n%s\n" % (n, "".join("ACGTN"[c] for c in tiny.ref.contig_codes(i))) for i, n in enumerate(tiny.ref.names))
	assert r.stdout == want
	for (name, seq), line in zip(tiny.contigs, r.stdout.split("\n")[1::2]):
		assert line == "".join(c.upper() if c in "ACGTacgt" else "N" for c in seq)


def test_getref_edge_and_usage(edge):
	"""getref on contigs of 0 and 1 bp; without arguments it prints its usage and fails."""
	contigs, prefix, ref = edge
	r = H.run_minibwa(["getref", prefix + ".l2b"], check=True)
	assert r.stdout == "".join(">%s\n%s\n" % (n, "".join(c.upper() if c in "ACGTacgt" else "N" for c in s)) for n, s in contigs)
	r = H.run_minibwa(["getref"])
	assert r.returncode == 1 and "Usage: minibwa getref" in r.stderr and r.stdout == ""


# ---------------------------------------------------------------------------
# bench
# ---------------------------------------------------------------------------

_M64 = (1 << 64) - 1


def _signed(v):
	return v - (1 << 64) if v >= 1 << 63 else v


def _bench_expect(ref, kind, n, intv=20):
	"""bench's points: positions k drawn with splitmix64 seeded with 11; 2a prints the C
	count of BWT rows [0, k), sa the SA value of row k, msa the XOR of the SA values of
	rows [k, min(k+intv, n)). Returns (printed values, checksum, rows that include the
	sentinel row 0)."""
	occ = ref.occ_table()
	x, cs, out, sentinel = 11, 1, [], set()
	for i in range(n):
		x, z = A.splitmix64(x)
		k = z % ref.n
		if kind == "2a":
			x, z = A.splitmix64(x)
			out.append(occ[k][1])
			continue
		rows = [k] if kind == "sa" else list(range(k, min(k + intv, ref.n)))
		v = 0
		for row in rows:
			v ^= _M64 if row == 0 else ref.sa[row - 1]
		if 0 in rows:
			sentinel.add(i)
		out.append(_signed(v))
		cs = ((cs * 0xBF58476D1CE4E5B9) & _M64) ^ v
	return out, cs, sentinel


@pytest.mark.parametrize("args,kind,n,intv", [
	(["-b", "sa", "-n", "3k"], "sa", 3000, 20),
	(["-b", "2a", "-n", "1500"], "2a", 1500, 20),
	(["-n", "400"], "2a", 400, 20),
	(["-b", "msa", "-n", "800", "-v", "7"], "msa", 800, 7),
	(["-b", "msa", "-n", "600", "-1"], "msa", 600, 20),
])
def test_bench_points(tiny, args, kind, n, intv):
	"""minibwa bench -p prints one value per point: for -b 2a the rank of C at a random row,
	for -b sa a random row's SA value, for -b msa the XOR of the SA values of -v rows
	(batched, or one by one with -1); -n takes k/m suffixes. For sa and msa the checksum on
	stderr is recomputed too (2a's checksum also uses the second rank, which bench draws
	independently and may pass in decreasing order)."""
	r = H.run_minibwa(["bench", "-p"] + args + [tiny.index + ".mbw"], check=True)
	got = [int(v) for v in r.stdout.split()]
	want, cs, sentinel = _bench_expect(tiny.ref, kind, n, intv)
	assert len(got) == n
	for i, (a, b) in enumerate(zip(got, want)):
		if i not in sentinel:
			assert a == b, "point %d" % i
	m = re.search(r"checksum = ([0-9a-f]+)", r.stderr)
	assert m
	if kind != "2a" and not sentinel:
		assert int(m.group(1), 16) == cs


def test_bench_usage(tiny):
	r = H.run_minibwa(["bench", "--help"])
	assert r.returncode == 0 and "Usage: minibwa bench" in r.stdout
	r = H.run_minibwa(["bench"])
	assert r.returncode == 1 and "Usage: minibwa bench" in r.stderr
