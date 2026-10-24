"""Reference implementations and fixtures for the area B tests (chaining and alignment).

Everything here is written from the documented definitions, not from
minibwa's code paths:

	aln_score        affine-gap score of a CIGAR; with two gap models a gap of
	                 length l costs min(q + l*e, q2 + l*e2) (ksw2's dual affine)
	segment_score    best-scoring segment of an alignment (the `ms:i` tag): the
	                 running score restarts at 0, a gap costs q+e (short reads)
	                 or q + e*log2(1+l) (long reads)
	mlen_blen        PAF columns 10-11: matching bases and alignment length,
	                 ambiguous bases excluded from both
	compute_ds       the `ds:Z` string: cs with the part of each indel that
	                 repeats the adjacent sequence in brackets
	eqx_cigar        an M-CIGAR rewritten with '=' and 'X' (N is always 'X')
	chain_best       brute-force chaining DP over anchors (the minimap2 score)
	mb_log2          minibwa's float32 log2 approximation, reproduced bit by bit
	make_bgenome     a genome with designed repeats for the hand-worked cases
"""

import re
import struct

import pytest

import minibwa_helpers as H

# ---------------------------------------------------------------------------
# Scores
# ---------------------------------------------------------------------------

DEFAULT_SC = dict(a=2, b=8, q=12, e=2, q2=23, e2=1, b_ambi=1)


def base_score(r, q, a=2, b=8, b_ambi=1):
	x, y = H.nt4(r), H.nt4(q)
	if x > 3 or y > 3:
		return -b_ambi
	return a if x == y else -b


def gap_cost(l, q, e, q2, e2):
	return min(q + l * e, q2 + l * e2)


def aligned_query(rec, read_fwd):
	"""The read bases a SAM record covers (hard clips removed), in reference orientation."""
	s = read_fwd if not rec.is_reverse else H.revcomp(read_fwd)
	ops = rec.cigar_ops
	h5 = ops[0][0] if ops and ops[0][1] == "H" else 0
	h3 = ops[-1][0] if len(ops) > 1 and ops[-1][1] == "H" else 0
	return s[h5:len(s) - h3]


def aln_score(ref, pos0, seq, cigar, a=2, b=8, q=12, e=2, q2=23, e2=1, b_ambi=1):
	s = 0
	for x in H._walk(ref, pos0, seq, cigar):
		if x[0] == "M":
			s += base_score(x[1], x[2], a, b, b_ambi)
		elif x[0] in "ID":
			s -= gap_cost(len(x[1]), q, e, q2, e2)
	return s


def _f32(x):
	return struct.unpack("f", struct.pack("f", x))[0]


def mb_log2(x):
	"""mbpriv.h mb_log2(): exponent bits plus a quadratic on the mantissa, in float32."""
	x = _f32(x)
	i = struct.unpack("I", struct.pack("f", x))[0]
	log_2 = _f32(((i >> 23) & 255) - 128)
	i = (i & ~(255 << 23)) + (127 << 23)
	z = struct.unpack("f", struct.pack("I", i & 0xffffffff))[0]
	t = _f32(_f32(_f32(_f32(-0.34484843) * z) + _f32(2.02466578)) * z)
	t = _f32(t - _f32(0.67487759))
	return _f32(log_2 + t)


def segment_score(ref, pos0, seq, cigar, long_mode=False, a=2, b=8, q=12, e=2, b_ambi=1):
	"""ms:i: the best segment of the alignment, the running score floored at 0, rounded."""
	s = mx = 0.0
	for x in H._walk(ref, pos0, seq, cigar):
		if x[0] == "M":
			s += base_score(x[1], x[2], a, b, b_ambi)
			if s < 0:
				s = 0.0
			else:
				mx = max(mx, s)
		elif x[0] in "ID":
			s -= q + e * mb_log2(1.0 + len(x[1])) if long_mode else q + e
			if s < 0:
				s = 0.0
	return int(mx + .499)


def mlen_blen(ref, pos0, seq, cigar):
	m = b = 0
	for x in H._walk(ref, pos0, seq, cigar):
		if x[0] == "M":
			if H.nt4(x[1]) < 4 and H.nt4(x[2]) < 4:
				b += 1
				m += x[1].upper() == x[2].upper()
		elif x[0] in "ID":
			b += sum(1 for c in x[1] if H.nt4(c) < 4)
	return m, b


def to_m_cigar(cigar):
	"""'3=1X4=2I' -> '8M2I'."""
	out = []
	for n, op in H.parse_cigar(cigar):
		op = "M" if op in "=X" else op
		if out and out[-1][1] == op:
			out[-1][0] += n
		else:
			out.append([n, op])
	return "".join("%d%s" % (n, op) for n, op in out)


def eqx_cigar(ref, pos0, seq, cigar):
	"""Rewrite M runs as '=' (both bases ACGT and equal) and 'X' (anything else)."""
	out = []
	r, q = pos0, 0

	def push(n, op):
		if out and out[-1][1] == op:
			out[-1][0] += n
		else:
			out.append([n, op])
	for n, op in H.parse_cigar(cigar):
		if op in "M=X":
			for k in range(n):
				x, y = H.nt4(ref[r + k]), H.nt4(seq[q + k])
				push(1, "=" if x == y and x < 4 else "X")
			r += n
			q += n
		elif op in "IS":
			if op == "I":
				push(n, "I")
			else:
				push(n, "S")
			q += n
		elif op in "DN":
			push(n, op)
			r += n
		else:
			push(n, op)
	return "".join("%d%s" % (n, op) for n, op in out)


def _ds_indel(s, y, n):
	"""Bracket the part of s[y:y+n] repeating its flanks: the longest prefix equal to the
	bases after the indel (within s) and the longest suffix equal to the bases before it."""
	code = [H.nt4(c) for c in s]
	lr = 0
	while lr < n and y - lr - 1 >= 0 and code[y + n - lr - 1] == code[y - lr - 1]:
		lr += 1
	ll = 0
	while ll < n and y + n + ll < len(s) and code[y + n + ll] == code[y + ll]:
		ll += 1
	ind = "".join("acgtn"[c] for c in code[y:y + n])
	if ll + lr >= n:
		return "[" + ind + "]"
	out = ""
	if ll > 0:
		out += "[" + ind[:ll] + "]"
	out += ind[ll:n - lr]
	if lr > 0:
		out += "[" + ind[n - lr:] + "]"
	return out


def compute_ds(ref, pos0, seq, cigar):
	"""ds:Z value (without the prefix) for an alignment without clipping in cigar: seq is
	the aligned query only. Context never reaches outside the aligned region."""
	ops = [(n, op) for n, op in H.parse_cigar(cigar) if op not in "SH"]
	tspan = sum(n for n, op in ops if op in "M=XDN")
	t = ref[pos0:pos0 + tspan]
	out, run = [], 0
	qo = to = 0
	for n, op in ops:
		if op in "M=X":
			for k in range(n):
				x, y = H.nt4(t[to + k]), H.nt4(seq[qo + k])
				if x == y:
					run += 1
				else:
					if run:
						out.append(":%d" % run)
						run = 0
					out.append("*%s%s" % ("acgtn"[x], "acgtn"[y]))
			qo += n
			to += n
		else:
			if run:
				out.append(":%d" % run)
				run = 0
			if op == "I":
				out.append("+" + _ds_indel(seq, qo, n))
				qo += n
			elif op == "D":
				out.append("-" + _ds_indel(t, to, n))
				to += n
	if run:
		out.append(":%d" % run)
	return "".join(out)


def merge_cs_runs(cs):
	"""Merge adjacent identical-run operations (':100:387' -> ':487'). With --eqx minibwa
	writes one run per CIGAR operation, and an 'X' run of N against N counts as identical in
	cs, so a run can be split in two; the merged string is the canonical form."""
	out = []
	for tok in re.findall(r":[0-9]+|[^:]+", cs):
		if tok[0] == ":" and out and out[-1][0] == ":":
			out[-1] = ":%d" % (int(out[-1][1:]) + int(tok[1:]))
		else:
			out.append(tok)
	return "".join(out)


def cigar_problems(cigar):
	"""Structural problems of a CIGAR: zero-length or repeated adjacent operations, clips
	inside, an indel next to a clip or at an end."""
	ops = H.parse_cigar(cigar)
	errs = []
	for i, (n, op) in enumerate(ops):
		if n == 0:
			errs.append("zero-length %s" % op)
		if i > 0 and ops[i - 1][1] == op:
			errs.append("repeated %s" % op)
	core = [op for n, op in ops]
	for clip in "HS":  # [H][S] core [S][H]
		if core and core[0] == clip:
			core = core[1:]
		if core and core[-1] == clip:
			core = core[:-1]
	if "S" in core or "H" in core:
		errs.append("clip inside the alignment")
	if core and (core[0] in "ID" or core[-1] in "ID"):
		errs.append("leading or trailing indel")
	return errs


def indels_left_aligned(ref, pos0, seq, cigar):
	"""Problems with indel placement: an indel between two matches that could move one base
	left with the same alignment (its last base equals the base before it)."""
	ops = H.parse_cigar(cigar)
	errs = []
	r, q = pos0, 0
	for i, (n, op) in enumerate(ops):
		inner = 0 < i < len(ops) - 1 and ops[i - 1][1] in "M=X" and ops[i + 1][1] in "M=X"
		if op in "M=X":
			r += n
			q += n
		elif op in "IS":
			if op == "I" and inner and H.nt4(seq[q - 1]) == H.nt4(seq[q + n - 1]):
				errs.append("insertion at query %d can move left" % q)
			q += n
		elif op == "D":
			if inner and H.nt4(ref[r - 1]) == H.nt4(ref[r + n - 1]):
				errs.append("deletion at %d can move left" % r)
			r += n
	return errs


# ---------------------------------------------------------------------------
# Chaining
# ---------------------------------------------------------------------------

CHAIN_PEN = _f32(_f32(0.8) * 0.01 * 19)  # chain_gap_scale * 0.01 * min_len, passed as float


class Anchor(object):
	__slots__ = ("sid", "qpos", "tpos", "len", "ctg", "rev")

	def __init__(self, ctg, rev, sid, qpos, tpos, length):
		self.ctg, self.rev, self.sid, self.qpos, self.tpos, self.len = ctg, rev, sid, qpos, tpos, length


def parse_dbg_anchors(stderr, qlens, ctg_ids):
	"""{qname: [Anchor]} from the `AC` lines of --dbg-anchor (qs on the forward read, ts
	the first base): converted to last-base coordinates on the anchor's strand."""
	out = {}
	for line in stderr.splitlines():
		f = line.split("\t")
		if f[0] != "AC":
			continue
		name, qs, strand, ctg, ts, length = f[1], int(f[2]), f[3], f[4], int(f[5]), int(f[6])
		rev = strand == "-"
		qpos = qlens[name] - 1 - qs if rev else qs + length - 1
		out.setdefault(name, []).append(Anchor(ctg, rev, ctg_ids[ctg] << 1 | rev, qpos, ts + length - 1, length))
	return out


def pair_score(aj, ai, max_dist, bw, pen=CHAIN_PEN):
	"""Score of appending anchor ai to a chain ending at aj, or None if they cannot chain."""
	dq, dr = ai.qpos - aj.qpos, ai.tpos - aj.tpos
	if ai.sid != aj.sid or dq <= 0 or dr <= 0:
		return None
	if dq > max_dist + ai.len or dr >= max_dist + ai.len:
		return None
	dd = abs(dr - dq)
	if dd > bw:
		return None
	sc = min(ai.len, dr, dq)
	if dd:
		sc -= int(_f32(_f32(pen * _f32(dd)) + _f32(0.5 * mb_log2(dd + 1))))
	return sc


def chain_best(anchors, max_dist, bw):
	"""(score, [anchors]) of the best chain by brute-force DP over all predecessor pairs."""
	a = sorted(anchors, key=lambda x: (x.sid, x.tpos, x.qpos))
	max_dist = max(max_dist, bw)
	f, p = [], []
	for i, ai in enumerate(a):
		best, bj = ai.len, -1
		for j in range(i - 1, -1, -1):
			sc = pair_score(a[j], ai, max_dist, bw)
			if sc is not None and f[j] + sc > best:
				best, bj = f[j] + sc, j
		f.append(best)
		p.append(bj)
	if not a:
		return 0, []
	top = max(range(len(a)), key=lambda i: (f[i], i))
	chain = []
	i = top
	while i >= 0:
		chain.append(a[i])
		i = p[i]
	return f[top], chain[::-1]


def chain_mlen_blen(chain):
	mlen = blen = chain[0].len
	for x, y in zip(chain, chain[1:]):
		tl, ql = y.tpos - x.tpos, y.qpos - x.qpos
		blen += max(tl, ql)
		mlen += y.len if tl > y.len and ql > y.len else min(tl, ql)
	return mlen, blen


def chain_coords(chain, qlen):
	"""(qs, qe, ts, te) of a chain, query coordinates on the forward read."""
	a0, a1 = chain[0], chain[-1]
	ts, te = a0.tpos + 1 - a0.len, a1.tpos + 1
	if not a0.rev:
		return a0.qpos + 1 - a0.len, a1.qpos + 1, ts, te
	return qlen - (a1.qpos + 1), qlen - (a0.qpos + 1 - a0.len), ts, te


def adap_value(x0, x1, length):
	"""The man page's adaptive value: x1 - (x1 - x0) * 2^(-(max(L,100)-100)/1900), rounded."""
	b = 2.0 ** (-(max(length, 100) - 100) / 1900.0)
	return int(x1 - (x1 - x0) * b + .499)


# ---------------------------------------------------------------------------
# Running and parsing
# ---------------------------------------------------------------------------


def body_records(res):
	return res.sam().records


def map_reads(index, path_dir, reads, args=(), name="b.fa", threads=None):
	"""Write (name, seq) tuples as FASTA and map them; returns the Result."""
	fa = H.write_fasta(str(path_dir / name), reads)
	return H.run_map(index, fa, list(args), threads=threads)


def one_record(res, qname=None):
	recs = [r for r in res.sam().records if qname is None or r.qname == qname]
	assert len(recs) == 1, "\n".join(r.line[:200] for r in recs)
	return recs[0]


def flip_base(c):
	"""A base different from c (deterministic)."""
	return {"A": "C", "C": "G", "G": "T", "T": "A"}[c.upper()]


# ---------------------------------------------------------------------------
# The designed genome
# ---------------------------------------------------------------------------


class BGenome(object):
	"""Contig b1 (about 51 kbp) and b2 (20 kbp) of random sequence with, on b1:
	a 1000 bp unit A at pos['A'] and a copy differing at offset 500 at pos['A2'];
	a 400 bp unit R three times (pos['R1'], pos['R2'] forward, pos['R3'] reverse-
	complemented) and once with a substitution at offset 200 (pos['R4']).
	b2 has a 50 bp N run at pos['N']."""

	def __init__(self, seed=8001):
		rng = H.Rng(seed)
		unit_a = rng.dna(1000)
		unit_a2 = unit_a[:500] + flip_base(unit_a[500]) + unit_a[501:]
		unit_r = rng.dna(400)
		unit_r4 = unit_r[:200] + flip_base(unit_r[200]) + unit_r[201:]
		parts = [("bg", rng.dna(5000)), ("A", unit_a), ("bg", rng.dna(5000)), ("A2", unit_a2), ("bg", rng.dna(5000)),
			("R1", unit_r), ("bg", rng.dna(4000)), ("R2", unit_r), ("bg", rng.dna(4000)), ("R3", H.revcomp(unit_r)),
			("bg", rng.dna(4000)), ("R4", unit_r4), ("bg", rng.dna(20000))]
		self.pos = {}
		seq, off = [], 0
		for name, s in parts:
			if name != "bg":
				self.pos[name] = off
			seq.append(s)
			off += len(s)
		b1 = "".join(seq)
		left, right = rng.dna(10000), rng.dna(9950)
		self.pos["N"] = len(left)
		b2 = left + "N" * 50 + right
		self.unit_a, self.unit_a2, self.unit_r = unit_a, unit_a2, unit_r
		self.genome = H.Genome([("b1", b1), ("b2", b2)])
		self.unique_start = self.pos["R4"] + 400 + 1000  # 18.6 kbp of plain random sequence follow

	def __getitem__(self, name):
		return self.genome[name]


@pytest.fixture(scope="module")
def bgen(index_of):
	g = BGenome()
	g.index = index_of(g.genome)
	return g
