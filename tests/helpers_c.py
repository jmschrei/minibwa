"""Helpers for area C (paired-end and Hi-C) tests.

Read pairs here are built with exact geometry instead of being simulated:
each mate is a stretch of a contig, optionally reverse-complemented, with
substitutions at chosen offsets. That gives exact 5' positions, so insert
sizes and pair orientations are known to the base, and lets a test plant
copies of a sequence (repeats, near-copies) where it needs them.

Coordinates follow pe.c's `mb_insert_dir`: the 5' position of a mate is its
leftmost base (0-based) on the forward strand and its exclusive end on the
reverse strand, and the distance of a pair is the absolute difference of the
two 5' positions. For an FR pair that is the fragment length.

Contents:
	Mate, Pair             a read with its true placement; a pair of them
	PlantedGenome          random contigs that tests can plant sequence into
	mate_at / pair_at      build mates and pairs from a genome
	fr_pair                an FR pair from a fragment start and length
	unseedable_subs        substitution offsets that leave no exact 19-mer
	write_pairs_c          read 1 and read 2 files (FASTQ)
	write_interleaved_c    one interleaved file, optionally with single reads
	reorient               turn a simulated FR pair into an FF, RR or RF pair, with its true code
	insert_dir             orientation code and 5'-to-5' distance (FF=0, FR=1, RF=2, RR=3)
	pestat                 the insert-size statistics minibwa infers, recomputed
	pair_score             the score of a pairing, recomputed
	stderr_pestat          parse the [M::mb_pestat] lines of a run's stderr
"""

import math
import os
import re

import minibwa_helpers as H

RLEN = 150
_SUB = {"A": "C", "C": "G", "G": "T", "T": "A", "a": "c", "c": "g", "g": "t", "t": "a"}


def substitute(seq, offsets):
	"""seq with the base at each offset replaced by a different base (A->C->G->T->A); other
	characters (N) are left alone."""
	s = list(seq)
	for o in offsets:
		s[o] = _SUB.get(s[o], s[o])
	return "".join(s)


def unseedable_subs(length=RLEN, step=15, first=7):
	"""Offsets of substitutions every `step` bases: the longest exact match is step-1 bases,
	shorter than minibwa's minimum seed length (19), so the read has no seed."""
	return list(range(first, length, step))


def _qual(name, n):
	"""A deterministic, varied quality string (so QUAL reversal is visible)."""
	h = 0
	for c in name:
		h = (h * 131 + ord(c)) & 0xffffffff
	rng = H.Rng(h)
	return "".join(chr(33 + 2 + rng.below(39)) for _ in range(n))


class Mate(object):
	"""One read and its true placement. seq/qual are as sequenced.

	contig, ts, te (0-based, exclusive) and rev give the alignment; subs are
	substitution offsets relative to ts in reference orientation. p5 is the
	5' position in pe.c's convention (te on the reverse strand)."""

	def __init__(self, name, seq, qual, contig, ts, te, rev, subs):
		self.name, self.seq, self.qual = name, seq, qual
		self.contig, self.ts, self.te, self.rev, self.subs = contig, ts, te, rev, list(subs)

	@property
	def p5(self):
		return self.te if self.rev else self.ts

	@property
	def strand(self):
		return "-" if self.rev else "+"

	@property
	def pos(self):
		return self.ts

	@property
	def fwd_seq(self):
		return H.revcomp(self.seq) if self.rev else self.seq

	@property
	def nm(self):
		return len(self.subs)

	def record(self):
		return (self.name, self.seq, self.qual)

	def __len__(self):
		return len(self.seq)

	def __repr__(self):
		return "Mate(%s %s:%d-%d%s subs=%d)" % (self.name, self.contig, self.ts, self.te, self.strand, len(self.subs))


class Pair(object):
	def __init__(self, name, r1, r2):
		self.name, self.r1, self.r2 = name, r1, r2

	def __iter__(self):
		return iter((self.r1, self.r2))

	@property
	def dir_dist(self):
		return insert_dir(self.r1.contig, self.r1.p5, self.r1.rev, self.r2.contig, self.r2.p5, self.r2.rev)

	def __repr__(self):
		return "Pair(%s, %r, %r)" % (self.name, self.r1, self.r2)


def mate_at(genome, name, contig, ts, rev, length=RLEN, subs=(), seq=None):
	"""A mate whose forward-strand sequence is genome[contig][ts:ts+length] with substitutions,
	or `seq` (forward strand) if given."""
	ref = genome[contig][ts:ts + length] if seq is None else seq
	if len(ref) != length:
		raise ValueError("mate %s at %s:%d runs off the contig" % (name, contig, ts))
	fwd = substitute(ref, subs)
	s = H.revcomp(fwd) if rev else fwd
	return Mate(name, s, _qual(name + ("/r" if rev else "/f"), length), contig, ts, ts + length, rev, subs)


def pair_at(genome, name, c1, p1, rev1, c2, p2, rev2, len1=RLEN, len2=RLEN, subs1=(), subs2=()):
	"""A pair from the 5' positions of its mates (pe.c convention: te for a reverse mate)."""
	m1 = mate_at(genome, name, c1, p1 - len1 if rev1 else p1, rev1, len1, subs1)
	m2 = mate_at(genome, name, c2, p2 - len2 if rev2 else p2, rev2, len2, subs2)
	return Pair(name, m1, m2)


def fr_pair(genome, name, contig, left, insert, r1_rev=False, len1=RLEN, len2=RLEN, subs1=(), subs2=()):
	"""An FR pair covering [left, left+insert): the forward mate starts at left, the reverse
	mate ends at left+insert. r1_rev=True makes read 1 the reverse mate."""
	if not r1_rev:
		return pair_at(genome, name, contig, left, False, contig, left + insert, True, len1, len2, subs1, subs2)
	return pair_at(genome, name, contig, left + insert, True, contig, left, False, len1, len2, subs1, subs2)


class PlantedGenome(object):
	"""Random contigs (splitmix64 from minibwa_helpers) that tests can write sequence into.

	plant(contig, pos, seq) overwrites; random(n) draws fresh sequence from a
	separate stream. genome() freezes the result as an H.Genome."""

	def __init__(self, seed, contigs):
		self.rng = H.Rng(seed)
		self.names = [n for n, _ in contigs]
		self.seqs = dict((n, list(self.rng.dna(l))) for n, l in contigs)
		self.extra = self.rng.fork("extra")
		self.planted = []

	def random(self, n):
		return self.extra.dna(n)

	def plant(self, contig, pos, seq):
		for c, p, l in self.planted:
			if c == contig and p < pos + len(seq) and pos < p + l:
				raise ValueError("plant at %s:%d overlaps an earlier plant at %d" % (contig, pos, p))
		self.seqs[contig][pos:pos + len(seq)] = list(seq)
		self.planted.append((contig, pos, len(seq)))

	def __getitem__(self, contig):
		return "".join(self.seqs[contig])

	def genome(self):
		return H.Genome([(n, "".join(self.seqs[n])) for n in self.names])


def write_pairs_c(d, pairs, tag="r", suffix=False):
	"""read 1 and read 2 FASTQ files for pairs; returns the two paths."""
	p1 = os.path.join(str(d), "%s_1.fq" % tag)
	p2 = os.path.join(str(d), "%s_2.fq" % tag)
	H.write_fastx(p1, [p.r1.record() for p in pairs], "fastq", suffix="/1" if suffix else None)
	H.write_fastx(p2, [p.r2.record() for p in pairs], "fastq", suffix="/2" if suffix else None)
	return p1, p2


def write_interleaved_c(path, items, suffix=False):
	"""One FASTQ file; items are Pairs (both mates, adjacent) or Mates (single reads)."""
	parts = []
	for it in items:
		if isinstance(it, Pair):
			parts.append(H.fastx_text([it.r1.record()], "fastq", "/1" if suffix else None))
			parts.append(H.fastx_text([it.r2.record()], "fastq", "/2" if suffix else None))
		else:
			parts.append(H.fastx_text([it.record()], "fastq"))
	with open(str(path), "w") as f:
		f.write("".join(parts))
	return str(path)


# Mates reverse-complemented to turn a simulated FR pair (read 1 at one end of the fragment, read 2
# at the other, facing each other) into another library orientation: FF flips read 2, RR read 1, RF both.
ORIENT_FLIPS = {"FR": (False, False), "FF": (False, True), "RR": (True, False), "RF": (True, True)}


def reorient(pair, orient):
	"""(read 1 seq, read 2 seq, (code, distance)) of an H.SimPair turned into orientation `orient`;
	the code and 5'-to-5' distance are those mb_pestat assigns from the true placement
	(insert_dir), on the read-1-first convention of pe.c."""
	out = []
	for r, flip in zip((pair.r1, pair.r2), ORIENT_FLIPS[orient]):
		rev = (r.strand == "-") != flip
		out.append((H.revcomp(r.seq) if flip else r.seq, rev, r.end if rev else r.pos))
	(s1, rev1, p1), (s2, rev2, p2) = out
	return s1, s2, insert_dir(pair.r1.contig, p1, rev1, pair.r2.contig, p2, rev2)


# ---------------------------------------------------------------------------
# Reference computations
# ---------------------------------------------------------------------------

DIR_NAMES = ("FF", "FR", "RF", "RR")


def insert_dir(c0, p0, rev0, c1, p1, rev1):
	"""(orientation code, 5'-to-5' distance) of read 1 at p0 and read 2 at p1 (5' positions).

	Codes as in pe.c: 1 (FR) when the leftmost mate is forward and the other
	reverse, 2 (RF) for the opposite, 0 (FF) when both mates point the same
	way and read 1 comes first in that direction, 3 (RR) otherwise. The
	contigs must match (None is returned otherwise)."""
	if c0 != c1:
		return None
	code = (int(rev0) << 1 | int(rev1)) ^ (0 if p0 < p1 else 3)
	return code, abs(p0 - p1)


def _int_trunc(x):
	return int(x)  # truncation toward zero, as a C cast


def pestat(dists_by_dir, min_dir_cnt=20, min_dir_ratio=0.05):
	"""Insert-size statistics per orientation from the distances of unique pairs.

	Returns four dicts (FF, FR, RF, RR) with failed, and when not failed
	p25/p50/p75, avg, std (over pairs within the outlier bounds), lo and hi.
	Follows the bwa-mem method that minibwa's mb_pestat keeps: quartiles by
	index round(q*n); outliers beyond 2 IQR excluded from the mean and
	standard deviation; proper-pair bounds 3 IQR beyond the quartiles,
	widened to at least 4 standard deviations around the mean, lo >= 1.
	"""
	ns = [len(dists_by_dir[d]) for d in range(4)]
	mx = max(ns)
	out = []
	for d in range(4):
		a = sorted(dists_by_dir[d])
		n = len(a)
		if n < min_dir_cnt or n < mx * min_dir_ratio:
			out.append({"failed": True, "n": n})
			continue
		p25 = a[_int_trunc(.25 * n + .499)]
		p50 = a[_int_trunc(.50 * n + .499)]
		p75 = a[_int_trunc(.75 * n + .499)]
		lo = _int_trunc(p25 - 2.0 * (p75 - p25) + .499)
		if lo < 1:
			lo = 1
		hi = _int_trunc(p75 + 2.0 * (p75 - p25) + .499)
		avg, x = 0.0, 0
		for v in a:
			if lo <= v <= hi:
				avg += v
				x += 1
		avg /= x
		std = 0.0
		for v in a:
			if lo <= v <= hi:
				std += (v - avg) * (v - avg)
		std = math.sqrt(std / x)
		lo = _int_trunc(p25 - 3.0 * (p75 - p25) + .499)
		hi = _int_trunc(p75 + 3.0 * (p75 - p25) + .499)
		if lo > avg - 4.0 * std:
			lo = _int_trunc(avg - 4.0 * std + .499)
		if hi < avg + 4.0 * std:
			hi = _int_trunc(avg + 4.0 * std + .499)
		if lo < 1:
			lo = 1
		n_lo = sum(1 for v in a if v < lo)
		n_hi = sum(1 for v in a if v > hi)
		out.append({"failed": False, "n": n, "p25": p25, "p50": p50, "p75": p75, "avg": avg, "std": std,
			"lo": lo, "hi": hi, "pct_lo": 100.0 * n_lo / n, "pct_hi": 100.0 * n_hi / n})
	return out


def log_term(dist, avg, std, match_sc=2):
	"""The insert-size term of a pair score: log of the two-sided normal tail probability
	of the distance, scaled by 1/log(4) and the match score."""
	ns = (dist - avg) / std
	t = 2.0 * math.erfc(abs(ns) * 0.707106781186547524401)
	if t <= 0.0:
		return float("-inf")
	return .721 * math.log(t) * match_sc


def pair_score(dp0, dp1, dist, avg, std, match_sc=2):
	"""Integer pair score: the two alignment scores plus log_term, clamped at 0, rounded."""
	s = dp0 + dp1 + log_term(dist, avg, std, match_sc)
	if s < 0.0:
		s = 0.0
	return int(s + .499)


# The insert-size report of mb_pestat on stderr is not documented, so it is parsed for its
# numbers only: a line naming the four orientations followed by four counts; per orientation, a
# line tagged ::FF/::FR/::RF/::RR with three integers then two decimals (quartiles; mean, std)
# and one with two integers then two decimals (bounds; percentages below/above). The words
# between the numbers are not matched, so rewording the messages does not break these tests.
_INT = r"(-?\d+)"
_DEC = r"(-?[0-9]+(?:\.[0-9]+)?)"
_PESTAT_CAND = re.compile(r"^\[M::[^\]]*\][^(\n]*\(FF, FR, RF, RR\)\D*" + r"\D+".join([_INT] * 4) + r"\D*$", re.M)
_PESTAT_PCT = re.compile(r"^\[M::[^\]\n]*::(FF|FR|RF|RR)\][^\n]*\(" + r", ".join([_INT] * 3) + r"\)[^(\n]*\(" + _DEC + r", " + _DEC +
	r"\)[^(\n]*$", re.M)
_PESTAT_BND = re.compile(r"^\[M::[^\]\n]*::(FF|FR|RF|RR)\][^(\n]*\(" + _INT + r", " + _INT + r"\)[^\d\n]*" + _DEC + r"/" + _DEC +
	r"[^\n]*$", re.M)


def stderr_pestat(stderr):
	"""Parsed [M::...] insert-size messages: a list (one per batch) of candidate counts (FF, FR, RF, RR),
	and lists of (dir, (p25, p50, p75, "avg", "std")) and (dir, (lo, hi, "pct_lo", "pct_hi")) over all
	batches, the decimals kept as printed."""
	cands = [tuple(int(x) for x in m.groups()) for m in _PESTAT_CAND.finditer(stderr)]
	pct = [(m.group(1), (int(m.group(2)), int(m.group(3)), int(m.group(4)), m.group(5), m.group(6))) for m in _PESTAT_PCT.finditer(stderr)]
	bnd = [(m.group(1), (int(m.group(2)), int(m.group(3)), m.group(4), m.group(5))) for m in _PESTAT_BND.finditer(stderr)]
	return cands, pct, bnd


def has_pestat_report(stderr):
	"""Whether stderr carries any insert-size statistics (a batch that inferred them)."""
	return stderr_pestat(stderr) != ([], [], [])
