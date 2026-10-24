"""Reference implementations for the area A tests (input and seeding).

Everything here is derived from the definitions, not from minibwa's code
paths: a suffix array of the indexed text built by prefix doubling, SA
intervals by binary search, occurrence counts over the BWT string, SMEMs and
sub-SMEMs by brute force, and the contig lookup of a text interval.

The indexed text. `minibwa index` encodes each FASTA base with A/C/G/T
(either case) as 0-3; every other character is ambiguous, recorded in the
.l2b ambiguity list and stored as a pseudo-random base drawn from splitmix64
seeded with 11 (one draw per ambiguous base, in file order). The BWT is built
over T = F + revcomp(F), where F is the concatenation of all contigs; SA row
0 is the sentinel, so the SA interval of a pattern P starts at
1 + #{suffixes of T lexicographically smaller than P}.
"""

import os

import minibwa_helpers as H

_M64 = (1 << 64) - 1
NT4 = dict((c, i) for i, c in enumerate("ACGT"))
NT4.update((c, i) for i, c in enumerate("acgt"))


def splitmix64(state):
	"""One step of splitmix64: (new state, output)."""
	state = (state + 0x9E3779B97F4A7C15) & _M64
	z = state
	z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & _M64
	z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & _M64
	return state, z ^ (z >> 31)


def nt4_codes(seq):
	"""minibwa's 2-bit-plus-N code of a read: A/C/G/T any case -> 0..3, anything else -> 4."""
	return [NT4.get(c, 4) for c in seq]


def revcomp_codes(codes):
	return [c if c > 3 else 3 - c for c in reversed(codes)]


def suffix_array(s):
	"""Suffix array of a str by prefix doubling (shorter suffix first on ties, as libsais)."""
	n = len(s)
	if n == 0:
		return []
	rank = [ord(c) for c in s]
	sa = list(range(n))
	k = 1
	while True:
		key = [(rank[i], rank[i + k] if i + k < n else -1) for i in range(n)]
		sa.sort(key=key.__getitem__)
		new = [0] * n
		r = 0
		for j in range(1, n):
			if key[sa[j]] != key[sa[j - 1]]:
				r += 1
			new[sa[j]] = r
		rank = new
		if r == n - 1:
			return sa
		k <<= 1


class RefText(object):
	"""The text minibwa indexes for a list of (name, sequence) contigs, with a naive suffix array.

	Attributes: names, lens, offs (contig offsets in F), tot (len(F)), fwd
	(F as an 'ACGT' str), text (T = F + revcomp(F)), n (len(T)), sa (suffix
	array of T, without the sentinel row), ambi (maximal runs [st, en) of
	ambiguous bases in F), amb (per-base ambiguity flags of F).
	"""

	def __init__(self, contigs, seed=11, with_sa=True):
		self.names = [c[0] for c in contigs]
		self.lens = [len(c[1]) for c in contigs]
		self.offs = []
		state = seed
		bases, amb = [], []
		off = 0
		for name, seq in contigs:
			self.offs.append(off)
			off += len(seq)
			for ch in seq:
				c = NT4.get(ch)
				if c is None:
					state, z = splitmix64(state)
					c = z & 3
					amb.append(1)
				else:
					amb.append(0)
				bases.append(c)
		self.tot = off
		self.amb = amb
		self.fwd = "".join("ACGT"[c] for c in bases)
		self.text = self.fwd + H.revcomp(self.fwd)
		self.n = len(self.text)
		self.ambi = []
		i = 0
		while i < self.tot:
			if amb[i]:
				j = i
				while j < self.tot and amb[j]:
					j += 1
				self.ambi.append((i, j))
				i = j
			else:
				i += 1
		self.sa = suffix_array(self.text) if with_sa else None
		self._occ = None

	# --- suffix-array intervals ---------------------------------------

	def _lower(self, p):
		t, sa, m = self.text, self.sa, len(p)
		lo, hi = 0, len(sa)
		while lo < hi:
			mid = (lo + hi) // 2
			if t[sa[mid]:sa[mid] + m] < p:
				lo = mid + 1
			else:
				hi = mid
		return lo

	def _upper(self, p, lo=0):
		t, sa, m = self.text, self.sa, len(p)
		hi = len(sa)
		while lo < hi:
			mid = (lo + hi) // 2
			if t[sa[mid]:sa[mid] + m] <= p:
				lo = mid + 1
			else:
				hi = mid
		return lo

	def interval(self, p):
		"""(x0, size) of pattern p (an 'ACGT' str) in minibwa's row numbering (row 0 = sentinel)."""
		lo = self._lower(p)
		return lo + 1, self._upper(p, lo) - lo

	def bi_interval(self, p):
		"""(x0, x1, size): the forward interval of p and of its reverse complement."""
		x0, size = self.interval(p)
		x1, size1 = self.interval(H.revcomp(p))
		assert size == size1
		return x0, x1, size

	def count(self, p):
		return self.interval(p)[1]

	def sa_row(self, row):
		"""SA value of a minibwa row (row 0 is the sentinel suffix at position n)."""
		return self.n if row == 0 else self.sa[row - 1]

	def positions(self, p):
		"""Text positions of all occurrences of p, in suffix-array order."""
		x0, size = self.interval(p)
		return [self.sa[x0 - 1 + j] for j in range(size)]

	# --- BWT ------------------------------------------------------------

	def bwt_full(self):
		"""The BWT with the sentinel: list of codes 0..3, None at the primary row."""
		out = []
		for row in range(self.n + 1):
			p = self.sa_row(row)
			out.append(None if p == 0 else "ACGT".index(self.text[p - 1]))
		return out

	def occ_table(self):
		"""occ[k] = (counts of A, C, G, T in BWT rows [0, k)) for k in 0..n+1."""
		if self._occ is None:
			c = [0, 0, 0, 0]
			occ = [tuple(c)]
			for b in self.bwt_full():
				if b is not None:
					c[b] += 1
				occ.append(tuple(c))
			self._occ = occ
		return self._occ

	def primary(self):
		return self.sa.index(0) + 1

	def l2(self):
		"""Cumulative base counts C(): L2[c] = number of bases < c in T."""
		cnt = [self.text.count(b) for b in "ACGT"]
		return [0, cnt[0], cnt[0] + cnt[1], cnt[0] + cnt[1] + cnt[2], self.n]

	# --- contig lookup ----------------------------------------------------

	def intv2cid(self, st, en):
		"""Contig of the text interval [st, en): (tid, cst, rev), or -1 (spans contigs),
		-2 (spans the strand boundary) or -3 (beyond the text)."""
		if en > 2 * self.tot:
			return -3
		if st < self.tot < en:
			return -2
		rev = 1 if st >= self.tot else 0
		s = st if not rev else 2 * self.tot - en
		for i, (o, l) in enumerate(zip(self.offs, self.lens)):
			if o <= s < o + l:
				return (i, s - o, rev) if s + (en - st) <= o + l else -1
		return -1

	def contig_codes(self, tid):
		"""l2b_getseq codes of a contig: 0..3, or 4 at ambiguous bases."""
		o, l = self.offs[tid], self.lens[tid]
		return [4 if self.amb[o + i] else "ACGT".index(self.fwd[o + i]) for i in range(l)]


# ---------------------------------------------------------------------------
# SMEMs by brute force
# ---------------------------------------------------------------------------


def codes_to_str(codes):
	return "".join("ACGTN"[c] for c in codes)


def right_extents(ref, q, st, en, min_occ):
	"""R[s] for s in [st, en): the largest e <= en with q[s:e] free of N and occurring at
	least min_occ times in the text (R[s] = s if no such e > s). Non-decreasing in s."""
	qs = codes_to_str(q)
	R = {}
	e = st
	for s in range(st, en):
		if e < s:
			e = s
		while e < en and q[e] < 4 and ref.count(qs[s:e + 1]) >= min_occ:
			e += 1
		R[s] = e
	return R


def smems(ref, q, min_len, min_occ=1, st=0, en=None):
	"""SMEMs of q (codes) within [st, en): maximal exact matches occurring at least min_occ
	times, at least min_len long, not contained in another such match starting in [st, en).
	Returns [(qs, qe, x0, x1, size)] ordered by qs."""
	if en is None:
		en = len(q)
	R = right_extents(ref, q, st, en, min_occ)
	qs = codes_to_str(q)
	out = []
	for s in range(st, en):
		e = R[s]
		if e - s >= min_len and (s == st or R[s - 1] < e):
			x0, x1, size = ref.bi_interval(qs[s:e])
			out.append((s, e, x0, x1, size))
	return out


def seed_intervals(ref, q, min_len, max_sub_occ):
	"""The intervals mb_seed_intv computes, as a set: the SMEMs (occurrence >= 1), then for
	each SMEM at least 2*min_len long occurring at most max_sub_occ times, the SMEMs inside
	it that occur more often than it, at least max(len/2, min_len) long."""
	out = set(smems(ref, q, min_len, 1))
	for (s, e, x0, x1, size) in list(out):
		if e - s < 2 * min_len or size > max_sub_occ:
			continue
		sub_len = max((e - s) // 2, min_len)
		out.update(smems(ref, q, sub_len, size + 1, s, e))
	return out


def sampled_rows(x0, size, max_occ):
	"""SA rows of an interval that become anchors: all of them, or max_occ rows spread
	evenly when the interval is larger than max_occ."""
	if size <= max_occ:
		return list(range(x0, x0 + size))
	rows, j, n = [], 0, 0
	while j < size and n < max_occ:
		step = max((size - j) // (max_occ - n), 1)
		rows.append(x0 + j)
		j += step
		n += 1
	return rows


def anchors(ref, qlen, seeds, max_occ):
	"""Anchors of a read from its seed intervals, as a set of (tid, rev, qs, ts, len): qs is
	the start on the read as given, ts the 0-based start on the contig's forward strand.
	Occurrences spanning two contigs or the strand boundary are dropped; an anchor contained
	in another on the same contig, strand and diagonal is dropped."""
	raw = set()
	for (s, e, x0, x1, size) in seeds:
		ln = e - s
		for row in sampled_rows(x0, size, max_occ):
			p = ref.sa_row(row)
			r = ref.intv2cid(p, p + ln)
			if not isinstance(r, tuple):
				continue
			tid, cst, rev = r
			raw.add((tid, rev, s, cst, ln))
	return drop_contained(raw, qlen)


def _diag(a, qlen):
	tid, rev, qs, ts, ln = a
	q = qlen - qs - ln if rev else qs  # start on the strand the anchor aligns to
	return ts - q


def drop_contained(raw, qlen):
	by = {}
	for a in raw:
		by.setdefault((a[0], a[1], _diag(a, qlen)), []).append(a)
	keep = set()
	for grp in by.values():
		for a in grp:
			contained = False
			for b in grp:
				if b != a and b[3] <= a[3] and a[3] + a[4] <= b[3] + b[4]:
					contained = True
					break
			if not contained:
				keep.add(a)
	return keep


# ---------------------------------------------------------------------------
# Test genomes and reads
# ---------------------------------------------------------------------------


def _iupac(rng, n):
	return "".join(rng.choice("RYKMSWBDHVN") for _ in range(n))


def make_tiny():
	"""A ~3.5 kbp genome for exhaustive checks against the naive suffix array.

	Contigs t1 (2000 bp, with a comment), t2 (1200 bp, starting with an N run, with a
	lowercase n two bases before its end), four short contigs of 7, 33, 64 and 129 bp and
	t3 (1500 bp). A 100 bp family F1 has three copies on both strands; a 40 bp family F2
	has 14 exact copies of its own; t3 holds ten copies of a 60 bp F4 that contains F2, so
	an F4 match occurs exactly 10 times (minibwa's sub-SMEM limit) and F2 24 times; there
	is a 300 bp lowercase stretch with an F1 copy in it, a 30 bp poly-A, a (CA)20 tract, a
	30 bp N run and scattered IUPAC codes. No contig ends in an ambiguous base: base minibwa leaves a
	contig's final run of ambiguous bases out of the .l2b ambiguity list, so l2b_getseq
	returns the random bases there instead of 4. Returns (contigs, features) with contigs
	a list of (name, seq) and features a dict of named positions."""
	rng = H.Rng(1701)
	t1 = list(rng.dna(2000))
	t2 = list(rng.dna(1200))
	f1 = rng.dna(100)
	f2 = rng.dna(40)
	feats = {"f1": f1, "f2": f2, "f1_pos": [], "f2_pos": []}

	def put(seq_list, pos, s, tag=None, ctg=None):
		seq_list[pos:pos + len(s)] = list(s)
		if tag:
			feats[tag].append((ctg, pos))

	put(t1, 100, f1, "f1_pos", "t1")
	put(t1, 700, H.revcomp(f1), "f1_pos", "t1")
	put(t1, 1100, "A" * 30)
	put(t1, 1200, "CA" * 20)
	put(t1, 1500, "N" * 30)
	for i, p in enumerate((260, 420, 560, 900, 1000, 1300, 1400, 1650, 1800, 1900)):
		put(t1, p, f2 if i % 3 else H.revcomp(f2), "f2_pos", "t1")
	put(t2, 50, f1.lower(), "f1_pos", "t2")
	for i, p in enumerate((200, 330, 600, 1000)):
		put(t2, p, f2 if i % 2 else H.revcomp(f2), "f2_pos", "t2")
	# lowercase stretch around the t2 copy of F1
	for i in range(20, 320):
		t2[i] = t2[i].lower()
	t2[0:5] = list("NNNNN")
	t2[-3] = "n"
	for p in (40, 333 + 50, 1700, 1960):
		t1[p] = _iupac(rng, 1)
	t2[800] = "R"
	t2[801] = "y"
	contigs = [("t1", "".join(t1)), ("t2", "".join(t2))]
	for name, n in (("s1", 7), ("s2", 33), ("s3", 64), ("s4", 129)):
		contigs.append((name, rng.dna(n)))
	# t3: ten copies of F4 = 10 bp + F2 + 10 bp, so an F4 match occurs exactly 10 times
	# (minibwa's sub-SMEM limit) and holds the more frequent F2
	f4 = rng.dna(10) + f2 + rng.dna(10)
	feats["f4"] = f4
	t3 = list(rng.dna(1500))
	for i in range(10):
		put(t3, 50 + 140 * i, f4 if i % 2 else H.revcomp(f4))
	contigs.append(("t3", "".join(t3)))
	return contigs, feats


def write_contigs(path, contigs, comments=None, line_width=60):
	comments = comments or {}
	with open(path, "w") as f:
		for name, seq in contigs:
			f.write(">%s%s\n" % (name, " " + comments[name] if name in comments else ""))
			for i in range(0, len(seq), line_width):
				f.write(seq[i:i + line_width] + "\n")
			if not seq:
				f.write("\n")
	return path


TINY_COMMENTS = {"t1": "first contig\twith a tab"}


def mutate_read(rng, s, sub_rate=0.02, indel=True):
	out = []
	i = 0
	while i < len(s):
		c = s[i]
		x = rng.random()
		if c.upper() in "ACGT" and x < sub_rate:
			out.append(rng.choice([b for b in "ACGT" if b != c.upper()]))
		elif indel and x < sub_rate * 1.15 and 10 < i < len(s) - 10:
			if rng.below(2):
				out.append(c)
				out.append(rng.dna(1))
			# else: deletion of c
		else:
			out.append(c)
		i += 1
	return "".join(out)


def tiny_reads(contigs, feats, ref):
	"""~180 reads that exercise seeding on the tiny genome: simulated reads from both
	strands, repeat copies, junctions between contigs and between the strands of the text,
	reads running past the start and the end of the text, N and IUPAC bases in reads,
	lowercase reads, reads through reference N runs (taken from the text, which has random
	bases there), homopolymer and dinucleotide reads, odd-length SMEMs with a sub-SMEM of
	exactly half their length, and reads shorter than a seed. Returns [(name, seq)]."""
	rng = H.Rng(4242)
	g = dict(contigs)
	out = []

	def add(name, seq):
		out.append(("%s_%d" % (name, len(out)), seq))

	for i in range(110):
		c = rng.choice(["t1", "t1", "t2", "s4", "s3"])
		L = rng.randint(30, 150)
		if len(g[c]) <= L:
			L = len(g[c]) - 1
		st = rng.below(len(g[c]) - L + 1)
		s = g[c][st:st + L]
		s = "".join(b if b.upper() in "ACGT" else rng.dna(1) for b in s)
		s = mutate_read(rng, s, sub_rate=rng.choice([0.0, 0.01, 0.03]))
		if rng.below(2):
			s = H.revcomp(s)
		add("sim", s)
	for k in range(6):
		u = mutate_read(rng, feats["f1"], 0.02, indel=False)
		add("f1", rng.dna(rng.randint(0, 30)) + (u if k % 2 else H.revcomp(u)) + rng.dna(rng.randint(0, 30)))
	for k in range(6):
		add("f2", rng.dna(rng.randint(10, 40)) + (feats["f2"] if k % 2 else H.revcomp(feats["f2"])) + rng.dna(rng.randint(10, 40)))
	# exact text substrings crossing contig ends and the strand boundary
	T, tot = ref.text, ref.tot
	for a, b in ((ref.offs[1], 40), (ref.offs[3], 25), (ref.offs[4], 30), (ref.offs[5], 50), (tot, 45), (tot, 12)):
		add("junc", T[a - b:a + b])
	# through the reference N run of t1 (the text has random bases there)
	add("nrun", T[1480:1560])
	add("nrun", H.revcomp(T[1490:1570]))
	add("nrun_asis", g["t1"][1470:1580])
	# N, IUPAC and lowercase in reads
	for k in range(6):
		s = list(T[rng.below(tot - 200):][:120])
		for _ in range(k + 1):
			s[rng.below(len(s))] = rng.choice("NNRYn.")
		add("ambig", "".join(s))
	add("lower", g["t1"][300:420].lower())
	add("mixed", g["t1"][600:700].lower() + g["t1"][700:760])
	add("polyA", "A" * 45)
	add("polyA_flank", g["t1"][1060:1100] + "A" * 30 + g["t1"][1130:1160])
	add("ca", "CA" * 30)
	add("ca_flank", g["t1"][1170:1200] + "CA" * 20 + g["t1"][1240:1280])
	add("short", g["t1"][500:510])
	add("short", g["t2"][500:518])
	add("exact_s1", g["s1"])
	add("exact_s2", H.revcomp(g["s2"]))
	add("chimera", g["t1"][200:260] + H.revcomp(g["t2"][900:960]))
	# an F4 match occurring exactly 10 times, bounded by the read ends or by N
	add("f4", feats["f4"])
	add("f4", H.revcomp(feats["f4"]))
	add("f4", rng.dna(25) + "N" + feats["f4"] + "N" + rng.dna(25))
	add("f4", rng.dna(15) + "N" + H.revcomp(feats["f4"]) + "N" + rng.dna(40))
	add("random", rng.dna(80))
	# odd-length (41 bp) SMEMs whose sub-SMEM is exactly floor(41/2) = 20 bp long: a unique
	# flank, then the first 19 bases of a forward copy of F2 (the base before F2 repeats too)
	for p in (420, 560, 1000, 1650, 1800):
		add("oddsub", g["t1"][p - 22:p + 19])
	# past either end of the text T = F + revcomp(F): a base, then the start of T (the start
	# of t1); the end of T (revcomp of t1's start), then a base. Base minibwa extends a size-1
	# interval at the sentinel's row by one base that is not there
	# (known-bugs/a-rank2a-sentinel-row.md, r1-rank2a-text-start.md), so the base next to the
	# text end is never that one: these reads check everything else about the text ends.
	bug = ref.bwt_full()[ref.primary() + 1]
	for c in range(4):
		if c == bug:
			continue
		add("tstart", rng.dna(20) + "ACGT"[c] + T[:60])
		add("tend", T[-60:] + "ACGT"[3 - c] + rng.dna(20))
	return out


def write_reads_fasta(path, reads):
	with open(path, "w") as f:
		for name, seq in reads:
			f.write(">%s\n%s\n" % (name, seq))
	return path


# ---------------------------------------------------------------------------
# Output parsers
# ---------------------------------------------------------------------------


def parse_fastmap(text):
	"""minibwa fastmap output -> [(name, len, [(qs, qe, size, positions or None)])]."""
	out = []
	cur = None
	for line in text.split("\n"):
		if not line:
			continue
		f = line.split("\t")
		if f[0] == "SQ":
			cur = (f[1], int(f[2]), [])
		elif f[0] == "EM":
			pos = None if f[4:] == ["*"] else f[4:]
			cur[2].append((int(f[1]), int(f[2]), int(f[3]), pos))
		elif f[0] == "//":
			out.append(cur)
			cur = None
		else:
			raise AssertionError("unexpected fastmap line %r" % line)
	assert cur is None, "fastmap output does not end with //"
	return out


def dbg_lines(stderr, tag):
	"""Tab-separated fields of the debug lines with the given tag (SD or AC) on stderr."""
	return [l.split("\t")[1:] for l in stderr.split("\n") if l.startswith(tag + "\t")]


def batch_lines(stderr):
	"""(bases, sequences) of each '[M::worker_pipeline...] mapped N bp in M sequences' line."""
	import re
	return [(int(m.group(1)), int(m.group(2))) for m in re.finditer(r"mapped (\d+) bp in (\d+) sequences", stderr)]



def make_edge():
	"""Contigs with awkward lengths for the .l2b lookups: an empty contig, a 1 bp contig,
	contigs of 31, 32, 33 and 64 bp around the 32-base words of the packed sequence, a
	leading N run, IUPAC codes and lowercase bases (no contig ends in an ambiguous base)."""
	rng = H.Rng(77)
	e4 = list(rng.dna(32))
	e4[5:9] = list("RYKM")
	e5 = list(rng.dna(33).lower())
	e5[16] = "n"
	return [("e1", "G"), ("e0", ""), ("e2", "NNNNNACGTACGTTTGACCA"), ("e3", rng.dna(31)), ("e4", "".join(e4)),
		("e5", "".join(e5)), ("e6", rng.dna(64)), ("e7", rng.dna(100))]


class Tiny(object):
	"""The tiny genome, its index, its naive suffix array and the seeding reads."""

	def __init__(self, index_of, directory):
		self.contigs, self.feats = make_tiny()
		self.fasta = write_contigs(os.path.join(directory, "tiny.fa"), self.contigs, TINY_COMMENTS)
		self.index = index_of(self.fasta)
		self.ref = RefText(self.contigs)
		self.reads = tiny_reads(self.contigs, self.feats, self.ref)
		self.reads_fa = write_reads_fasta(os.path.join(directory, "tiny_reads.fa"), self.reads)
		self.codes = dict((n, nt4_codes(s)) for n, s in self.reads)


def parse_rows(text, tag):
	"""Driver output lines starting with tag, split on tabs (tag dropped)."""
	out = []
	for line in text.split("\n"):
		if line.startswith(tag + "\t"):
			out.append(line.split("\t")[1:])
	return out


def intervals_by_read(text, tag):
	"""{read name: set of (qs, qe, x0, x1, size)} from M/m/D lines of drv_a_bwt."""
	out = {}
	for f in parse_rows(text, tag):
		out.setdefault(f[0], set()).add(tuple(int(x) for x in f[1:]))
	return out


def warnings(stderr):
	"""stderr lines other than minibwa's [M::...] progress lines. Tests check that a warning
	appears and what it names, not its wording, which no minibwa document fixes."""
	return [l for l in stderr.splitlines() if l.strip() and not l.startswith("[M::")]
