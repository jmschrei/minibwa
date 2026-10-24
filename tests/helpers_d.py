"""Inputs and reference computations for the area D tests (output, CLI, API).

Everything here is deterministic (driven by minibwa_helpers.Rng) and uses
the standard library only.

	graded_repeat_genome()  two contigs carrying six copies of one 400 bp unit,
	                        each copy differing from the test read at a known
	                        number of positions, so every secondary hit has a
	                        score, NM and position known in advance
	mixed_reads(genome)     unique, repeat, chimeric, unmapped, short, long,
	                        IUPAC/lowercase and N-run reads with comments
	bisulfite_reads/_pairs  directional BS-seq reads with recorded truth
	ungapped_scores()       AS and ms of an ungapped alignment under minibwa's
	                        scoring (match 2, mismatch 8, ambiguous 1)
	expected_sa()           SA:Z entries recomputed from the other records
	mb_version()            MB_VERSION read from minibwa.h
"""

import os
import re

import minibwa_build
import minibwa_helpers as H

MATCH, MISMATCH, AMBIG = 2, 8, 1  # mb_opt_init(): a, b, b_ambi


def mb_version(src=None):
	"""MB_VERSION as written in minibwa.h of the source tree under test."""
	src = src or minibwa_build.src_dir()
	with open(os.path.join(src, "minibwa.h")) as f:
		m = re.search(r'^#define\s+MB_VERSION\s+"([^"]*)"', f.read(), re.M)
	if not m:
		raise AssertionError("no MB_VERSION in minibwa.h")
	return m.group(1)


def _sub(c):
	return "ACGT"[("ACGT".index(c) + 1) % 4]


# ---------------------------------------------------------------------------
# A genome whose secondary hits are known in advance
# ---------------------------------------------------------------------------

class GradedRepeats(object):
	"""Two contigs with two repeat families whose hits are known in advance.

	Family A: five copies of a 400 bp unit and a read from unit[100:250]. The
	read differs from every copy at read offsets 49 and 100 (so the exact
	match [50, 100) is a maximal seed in all five copies). Copy k also
	carries EXTRAS[k] substitutions at offsets from EXTRA_OFFSETS, all at
	least six bases apart and at least six from the read ends, so each copy
	aligns end to end without gaps. Under the default scoring a copy with m
	differences scores AS = ms = 300 - 10*m: 280, 260, 250, 210 and 200.

	Family B: seven identical copies of a 300 bp unit on both strands and a
	read from unit[50:200], so a read has six equally good secondary hits.
	"""

	EXTRAS = (0, 2, 3, 7, 8)
	EXTRA_OFFSETS = (106, 112, 118, 124, 130, 136, 142, 6)
	READ_DIFFS = (49, 100)
	PLAN_A = (("ctgA", 2000, "+"), ("ctgA", 9000, "-"), ("ctgB", 3000, "+"), ("ctgA", 16000, "+"), ("ctgB", 11000, "-"))
	PLAN_B = (("ctgA", 5000, "+"), ("ctgA", 12500, "-"), ("ctgA", 20000, "+"), ("ctgB", 6500, "-"),
		("ctgB", 14500, "+"), ("ctgB", 18000, "+"), ("ctgB", 20500, "-"))

	def __init__(self):
		rng = H.Rng(101)
		unit = rng.dna(400)
		unit_b = rng.dna(300)
		seqs = {"ctgA": list(rng.dna(24000)), "ctgB": list(rng.dna(22000))}
		read = list(unit[100:250])
		for o in self.READ_DIFFS:
			read[o] = _sub(read[o])
		self.read = "".join(read)
		self.read_b = unit_b[50:200]
		self.copies, self.copies_b = [], []
		feats = []
		for k, (ctg, p, strand) in enumerate(self.PLAN_A):
			cp = list(unit)
			for o in self.EXTRA_OFFSETS[:self.EXTRAS[k]]:
				cp[100 + o] = _sub(cp[100 + o])
			win = "".join(cp[100:250])
			cp = "".join(cp)
			if strand == "-":
				cp = H.revcomp(cp)
			seqs[ctg][p:p + 400] = list(cp)
			feats.append(H.Feature(ctg, p, p + 400, "repeat", (0, k, strand)))
			pos1 = p + 101 if strand == "+" else p + 151  # 1-based leftmost base of the read's alignment
			nm = sum(a != b for a, b in zip(win, self.read))
			self.copies.append({"ctg": ctg, "pos": pos1, "strand": strand, "nm": nm,
				"score": MATCH * (150 - nm) - MISMATCH * nm})
		for k, (ctg, p, strand) in enumerate(self.PLAN_B):
			seqs[ctg][p:p + 300] = list(unit_b if strand == "+" else H.revcomp(unit_b))
			feats.append(H.Feature(ctg, p, p + 300, "repeat", (1, k, strand)))
			self.copies_b.append({"ctg": ctg, "pos": p + 51 if strand == "+" else p + 101, "strand": strand, "nm": 0})
		feats.sort(key=lambda f: (f.contig, f.start))
		self.genome = H.Genome([(c, "".join(seqs[c])) for c in sorted(seqs)], feats)

	def reads(self):
		"""Family A's read and its reverse complement, then family B's read, with per-base qualities."""
		q = "".join(chr(35 + (i * 7) % 40) for i in range(150))
		return [("q", self.read, q), ("qr", H.revcomp(self.read), q[::-1]), ("b", self.read_b, q)]

	def secondaries_at(self, out_s):
		"""Family A copies 1..4 whose score passes `dp_max >= out_s * primary` (out_s as a C float)."""
		import struct
		f = struct.unpack("f", struct.pack("f", out_s))[0]
		best = self.copies[0]["score"]
		return [c for c in self.copies[1:] if c["score"] >= f * best]


# ---------------------------------------------------------------------------
# Read sets
# ---------------------------------------------------------------------------

def _qual(rng, n):
	return "".join(chr(33 + 2 + rng.below(39)) for _ in range(n))


def _unique_start(genome, rng, ctg, length, avoid=()):
	for _ in range(10000):
		p = rng.randint(300, genome.length(ctg) - length - 300)
		if genome.is_unique(ctg, p, p + length, 50) and all(abs(p - a) > 2000 for a in avoid):
			return p
	raise ValueError("no unique start")


def mixed_reads(genome, seed=41):
	"""A read set that exercises every branch of the SAM/PAF writers.

	Returns (records, info). records are (name, seq, qual, comment) tuples;
	info maps a read name to a dict with its category ('unique', 'repeat',
	'repeat_edge', 'chimeric', 'random', 'short', 'long', 'iupac', 'nrun') and, where known,
	its truth (contig, pos (0-based), strand; for chimeras a list 'parts'
	of (qstart, qend, contig, pos, strand)).
	"""
	rng = H.Rng(seed)
	recs, info = [], {}

	def add(name, seq, cat, comment=None, **truth):
		recs.append((name, seq, _qual(rng, len(seq)), comment))
		d = {"cat": cat}
		d.update(truth)
		info[name] = d

	sim = H.simulate_reads(genome, 40, length=150, seed=seed + 1, sub_rate=0.01, indel_rate=0.002, n_rate=0.005, unique=True,
		unique_pad=50)
	for i, r in enumerate(sim):
		add("u%d" % i, r.seq, "unique", "BC:Z:%s" % H.Rng(i).dna(8) if i % 3 == 0 else None,
			contig=r.contig, pos=r.pos, strand=r.strand)
	reps = [f for f in genome.features if f.kind == "repeat"]
	for i in range(8):
		f = reps[i % len(reps)]
		st = f.start + 20 + rng.below(f.end - f.start - 190)
		s = genome[f.contig][st:st + 150]
		add("rp%d" % i, s if i % 2 == 0 else H.revcomp(s), "repeat")
	for i in range(2):  # 100 bp of a repeat copy and 50 bp of unique flank: the other copies give clipped secondaries
		f = reps[i]
		s = genome[f.contig][f.end - 100:f.end + 50]
		add("re%d" % i, s if i % 2 == 0 else H.revcomp(s), "repeat_edge", contig=f.contig, pos=f.end - 100, strand="+-"[i % 2])
	for i in range(8):
		la = 60 + rng.below(31)
		lb = 150 - la
		ca = genome.names[i % 2]
		cb = genome.names[(i + 1) % 2] if i % 3 else ca
		pa = _unique_start(genome, rng, ca, la)
		pb = _unique_start(genome, rng, cb, lb, avoid=(pa,) if cb == ca else ())
		sa = "+" if i % 2 == 0 else "-"
		sb = "-" if i % 4 < 2 else "+"
		a = genome.fetch(ca, pa, pa + la, sa)
		b = genome.fetch(cb, pb, pb + lb, sb)
		add("ch%d" % i, a + b, "chimeric", "XX:i:%d" % i if i % 2 else None,
			parts=[(0, la, ca, pa, sa), (la, 150, cb, pb, sb)])
	for i in range(3):
		add("rnd%d" % i, rng.dna(120), "random")
	for i, r in enumerate(H.simulate_reads(genome, 3, length=45, seed=seed + 2, sub_rate=0.0, indel_rate=0.0, unique=True,
			unique_pad=50)):
		add("sh%d" % i, r.seq, "short", contig=r.contig, pos=r.pos, strand=r.strand)
	for i, r in enumerate(H.simulate_long_reads(genome, 2, length=(1500, 2500), seed=seed + 3, contigs=["chr1"])):
		add("lr%d" % i, r.seq, "long", contig=r.contig, pos=r.pos, strand=r.strand)
	iupac = "RYKMSWBDHVN"
	for i, r in enumerate(H.simulate_reads(genome, 4, length=150, seed=seed + 4, sub_rate=0.0, indel_rate=0.0, unique=True,
			unique_pad=50)):
		s = list(r.seq)
		for j in range(3):
			s[30 + 40 * j] = iupac[(i * 3 + j) % len(iupac)]
		s = "".join(s)
		if i % 2 == 0:
			s = s[:60].lower() + s[60:]
		add("iu%d" % i, s, "iupac", contig=r.contig, pos=r.pos, strand=r.strand)
	for i, f in enumerate([f for f in genome.features if f.kind == "N"]):
		st = f.start - 100
		s = genome[f.contig][st:st + 150]
		add("nr%d" % i, s if i % 2 == 0 else H.revcomp(s), "nrun", contig=f.contig, pos=st, strand="+-"[i % 2])
	return recs, info


def write_records(path, recs):
	"""FASTQ with comments; returns path."""
	return H.write_fastq(path, [(n, s, q, c) for n, s, q, c in recs])


# ---------------------------------------------------------------------------
# Bisulfite reads
# ---------------------------------------------------------------------------

def _convert(seq, rng, keep_cpg):
	"""C->T on one strand, keeping a C before G with probability keep_cpg (methylated CpG)."""
	out = list(seq)
	for i, c in enumerate(seq):
		if c == "C":
			if i + 1 < len(seq) and seq[i + 1] == "G" and rng.chance(keep_cpg):
				continue
			out[i] = "T"
	return "".join(out)


def bisulfite_pairs(genome, n, seed=61, read_len=150, insert=(250, 400), keep_cpg=0.7, sub_rate=0.0):
	"""Directional BS-seq pairs from unique sequence.

	For an original-top (OT) fragment the top strand is converted C->T, read 1
	is its first read_len bases and read 2 the reverse complement of its last
	read_len bases; for an original-bottom (OB) fragment the same is done on the
	bottom strand. In SAM orientation, OT reads show C->T against the
	reference and OB reads G->A. sub_rate adds substitutions after conversion.
	Returns [(name, r1, r2, info)] with info = contig, frag_start, frag_end,
	strand ('OT'/'OB').
	"""
	rng = H.Rng(seed)
	out = []
	names = [c for c in genome.names if genome.length(c) > 5000]
	for i in range(n):
		for _ in range(10000):
			ins = rng.randint(insert[0], insert[1])
			ctg = names[rng.below(len(names))]
			st = rng.randint(200, genome.length(ctg) - ins - 200)
			if genome.is_unique(ctg, st, st + ins, 20):
				break
		frag = genome[ctg][st:st + ins].upper()
		kind = "OT" if i % 2 == 0 else "OB"
		strand = frag if kind == "OT" else H.revcomp(frag)
		conv = _convert(strand, rng, keep_cpg)
		if sub_rate > 0:
			conv = H.mutate(rng, conv, sub_rate)
		r1 = conv[:read_len]
		r2 = H.revcomp(conv[-read_len:])
		out.append(("bs%d" % i, r1, r2, {"contig": ctg, "frag_start": st, "frag_end": st + ins, "strand": kind}))
	return out


# ---------------------------------------------------------------------------
# Reference computations
# ---------------------------------------------------------------------------

def _mat(rb, qb, meth=None):
	rb, qb = rb.upper(), qb.upper()
	if rb not in "ACGT" or qb not in "ACGT":
		return -AMBIG
	if rb == qb or (meth == "C2T" and rb == "C" and qb == "T") or (meth == "G2A" and rb == "G" and qb == "A"):
		return MATCH
	return -MISMATCH


def ungapped_scores(ref, pos0, seq, meth=None):
	"""(AS, ms) of an ungapped alignment of seq at ref[pos0:]: AS is the sum of the
	column scores, ms the best-scoring segment (running sum floored at 0), as
	mb_update_extra() computes dp_max0 for a CIGAR of one M operation."""
	s = best = tot = 0
	for k, qb in enumerate(seq):
		x = _mat(ref[pos0 + k], qb, meth)
		tot += x
		s += x
		if s < 0:
			s = 0
		elif s > best:
			best = s
	return tot, best


def clip_lengths(rec):
	"""(left, right) clip lengths (soft or hard) in reference orientation."""
	return rec.left_clip, rec.right_clip


def query_interval(rec, qlen):
	"""[qs, qe) of the aligned part on the read as sequenced (5'->3')."""
	left, right = clip_lengths(rec)
	if rec.is_reverse:
		return right, qlen - left
	return left, qlen - right


def expected_sa(records, rec, qlen):
	"""SA:Z entries minibwa should write on rec: one per other non-secondary record of the same read,
	as (rname, pos, strand, cigar, mapq, nm); CIGAR is 5' clip S, then min(query, ref span) M, the
	excess as I or D, then 3' clip S, in reference orientation."""
	out = []
	for o in records:
		if o is rec or o.is_secondary:
			continue
		qs, qe = query_interval(o, qlen)
		qspan = qe - qs
		tspan = o.ref_span
		m = min(qspan, tspan)
		clip5 = qlen - qe if o.is_reverse else qs  # leftmost clip in reference orientation
		clip3 = qs if o.is_reverse else qlen - qe
		cig = ""
		if clip5:
			cig += "%dS" % clip5
		if m:
			cig += "%dM" % m
		if qspan > tspan:
			cig += "%dI" % (qspan - tspan)
		elif tspan > qspan:
			cig += "%dD" % (tspan - qspan)
		if clip3:
			cig += "%dS" % clip3
		out.append((o.rname, o.pos, o.strand, cig, o.mapq, o.tags["NM"]))
	return out


def parse_sa(value):
	out = []
	for e in value.split(";"):
		if not e:
			continue
		f = e.split(",")
		out.append((f[0], int(f[1]), f[2], f[3], int(f[4]), int(f[5])))
	return out


def parse_xa(value):
	"""XA:Z value -> [(rname, strand, pos, cigar, nm)]."""
	out = []
	for e in value.split(";"):
		if not e:
			continue
		f = e.split(",")
		out.append((f[0], f[1][0], int(f[1][1:]), f[2], int(f[3])))
	return out
