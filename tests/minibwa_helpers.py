"""Shared helpers for the minibwa test suite.

Standard library only. Everything random here is driven by `Rng`, a
splitmix64 generator written with integer arithmetic, so a given seed gives
the same genome and the same reads on every Python version and platform.

Contents:
	Rng                         deterministic random numbers
	make_genome / Genome        synthetic genomes with repeats, N runs,
	                            low-complexity sequence and short contigs
	simulate_reads / _pairs /   reads with recorded truth (contig, position,
	_hic_pairs / _long_reads    strand, CIGAR, edits)
	write_fasta / write_fastq / FASTA/FASTQ writers, plain or gzip
	write_pairs / write_interleaved
	read_fastx / read_fasta     FASTA/FASTQ readers, plain or gzip
	index_genome                run `minibwa index`, cached per genome per session
	run_minibwa / run_map       run the binary with a timeout and no core dumps
	parse_sam / parse_paf       strict SAM and PAF parsers with typed tags
	compute_nm / compute_md /   recompute alignment tags from the genome
	compute_cs / cigar_*
	sam_invariant_errors        generic SAM consistency checks
	check_pinned / strip_pg     compare output against a stored pin
	run_c_case / assert_c_case  run a case of a C driver in tests/c/
"""

import difflib
import gzip
import hashlib
import io
import os
import re
import subprocess

import minibwa_build

TESTS_DIR = minibwa_build.TESTS_DIR
REPO_ROOT = minibwa_build.REPO_ROOT
DATA_DIR = os.path.join(REPO_ROOT, "test")
PINS_DIR = os.path.join(TESTS_DIR, "pins")

# ---------------------------------------------------------------------------
# Environment
# ---------------------------------------------------------------------------


def test_threads():
	"""Worker threads (-t) for tests that do not pin it: $MINIBWA_TEST_THREADS, default 1."""
	try:
		n = int(os.environ.get("MINIBWA_TEST_THREADS", "1"))
	except ValueError:
		raise ValueError("MINIBWA_TEST_THREADS must be an integer")
	if n < 1:
		raise ValueError("MINIBWA_TEST_THREADS must be at least 1")
	return n


def run_timeout():
	"""Timeout in seconds for one minibwa process: $MINIBWA_TEST_TIMEOUT, default 120."""
	return float(os.environ.get("MINIBWA_TEST_TIMEOUT", "120"))


# ---------------------------------------------------------------------------
# Deterministic random numbers
# ---------------------------------------------------------------------------

_M64 = (1 << 64) - 1


class Rng(object):
	"""splitmix64. Integer arithmetic only, so streams match across Python versions."""

	def __init__(self, seed):
		self.state = int(seed) & _M64

	def u64(self):
		self.state = (self.state + 0x9E3779B97F4A7C15) & _M64
		z = self.state
		z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & _M64
		z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & _M64
		return z ^ (z >> 31)

	def below(self, n):
		"""Uniform integer in [0, n)."""
		n = int(n)
		if n <= 0:
			raise ValueError("below() needs n > 0")
		lim = (1 << 64) - ((1 << 64) % n)
		while True:
			x = self.u64()
			if x < lim:
				return x % n

	def randint(self, lo, hi):
		"""Uniform integer in [lo, hi], both inclusive."""
		return lo + self.below(hi - lo + 1)

	def random(self):
		"""Uniform float in [0, 1) with 53 random bits; exact, hence portable."""
		return (self.u64() >> 11) * (1.0 / (1 << 53))

	def chance(self, p):
		return p > 0 and self.random() < p

	def choice(self, seq):
		return seq[self.below(len(seq))]

	def weighted_index(self, weights):
		"""Index i with probability weights[i] / sum(weights); integer weights."""
		tot = sum(weights)
		x = self.below(tot)
		for i, w in enumerate(weights):
			if x < w:
				return i
			x -= w
		raise AssertionError("unreachable")

	def normal_int(self, mean, sd):
		"""Approximately normal integer (Irwin-Hall sum of 12 uniforms)."""
		s = 0
		for _ in range(12):
			s += self.below(1 << 16)
		s -= 6 << 16
		return int(mean + int(s * sd / float(1 << 16)))

	def dna(self, n, alphabet="ACGT"):
		out = []
		while len(out) < n:
			x = self.u64()
			for _ in range(32):
				out.append(alphabet[x & 3])
				x >>= 2
		return "".join(out[:n])

	def fork(self, tag):
		"""An independent generator derived from this one's state and a tag."""
		h = hashlib.sha256(("%d:%s" % (self.state, tag)).encode()).digest()
		return Rng(int.from_bytes(h[:8], "little"))


# ---------------------------------------------------------------------------
# Sequences
# ---------------------------------------------------------------------------

_COMP = {}
for _a, _b in zip("ACGTUNRYKMSWBDHVacgtunrykmswbdhv", "TGCAANYRMKSWVHDBtgcaanyrmkswvhdb"):
	_COMP[_a] = _b
_COMP_TABLE = str.maketrans(_COMP)


def revcomp(seq):
	return seq.translate(_COMP_TABLE)[::-1]


def nt4(c):
	"""minibwa's base code: A/C/G/T (any case) -> 0..3, anything else -> 4."""
	return "ACGT".find(c.upper()) if c.upper() in "ACGT" else 4


def mutate(rng, seq, sub_rate):
	"""Substitutions only, at rate sub_rate; returns the new sequence."""
	out = list(seq)
	for i, c in enumerate(out):
		if c.upper() in "ACGT" and rng.chance(sub_rate):
			out[i] = _other_base(rng, c)
	return "".join(out)


def _other_base(rng, c):
	alts = [b for b in "ACGT" if b != c.upper()]
	return alts[rng.below(3)]


# ---------------------------------------------------------------------------
# Synthetic genomes
# ---------------------------------------------------------------------------


class Feature(object):
	"""An annotated interval [start, end) on a contig; kind is 'repeat', 'N' or 'lowcomp'."""
	__slots__ = ("contig", "start", "end", "kind", "info")

	def __init__(self, contig, start, end, kind, info=None):
		self.contig, self.start, self.end, self.kind, self.info = contig, start, end, kind, info

	def __repr__(self):
		return "Feature(%r, %d, %d, %r, %r)" % (self.contig, self.start, self.end, self.kind, self.info)


class Genome(object):
	"""Contigs in order plus annotated features."""

	def __init__(self, contigs, features=()):
		self.names = [n for n, _ in contigs]
		if len(set(self.names)) != len(self.names):
			raise ValueError("duplicate contig names")
		self.seqs = dict(contigs)
		self.features = list(features)

	def __getitem__(self, name):
		return self.seqs[name]

	def __contains__(self, name):
		return name in self.seqs

	def length(self, name):
		return len(self.seqs[name])

	@property
	def total_length(self):
		return sum(len(s) for s in self.seqs.values())

	def fetch(self, contig, start, end, strand="+"):
		s = self.seqs[contig][start:end]
		return s if strand == "+" else revcomp(s)

	def fasta_text(self, line_width=60):
		out = []
		for n in self.names:
			out.append(">%s\n" % n)
			s = self.seqs[n]
			if line_width and line_width > 0:
				for i in range(0, len(s), line_width):
					out.append(s[i:i + line_width] + "\n")
			else:
				out.append(s + "\n")
		return "".join(out)

	def digest(self):
		return hashlib.sha256(self.fasta_text().encode()).hexdigest()

	def write_fasta(self, path, line_width=60, gz=None):
		_write_text(path, self.fasta_text(line_width), gz)
		return path

	def features_overlapping(self, contig, start, end, kinds=None):
		return [f for f in self.features if f.contig == contig and f.start < end and start < f.end and
			(kinds is None or f.kind in kinds)]

	def is_unique(self, contig, start, end, pad=0):
		"""True if [start-pad, end+pad) overlaps no repeat, N run or low-complexity feature."""
		return not self.features_overlapping(contig, start - pad, end + pad)


def make_genome(seed=1, contig_lengths=(60000, 40000), names=None, repeat_families=0, repeat_len=(300, 1000),
		repeat_copies=3, repeat_divergence=0.0, soft_mask_repeats=False, n_runs=0, n_run_len=(10, 500),
		low_complexity=0, low_complexity_len=(50, 300), short_contigs=0, short_contig_len=(20, 200), margin=200):
	"""A deterministic random genome.

	Contigs are uniform random DNA. Then, in this order and without
	overlapping each other or coming within `margin` bp of a contig end:
	`repeat_families` families of `repeat_copies` copies each (a random unit of
	length drawn from `repeat_len`, each copy on a random strand with
	substitutions at `repeat_divergence`, lowercase if `soft_mask_repeats`);
	`low_complexity` tracts (homopolymer, di- or trinucleotide repeats);
	`n_runs` runs of N. Finally `short_contigs` extra contigs with lengths
	drawn from `short_contig_len` are appended. Every placed interval is
	recorded in Genome.features.
	"""
	rng = Rng(seed)
	if names is None:
		names = ["chr%d" % (i + 1) for i in range(len(contig_lengths))]
	if len(names) != len(contig_lengths):
		raise ValueError("names and contig_lengths differ in length")
	seqs = [list(rng.dna(n)) for n in contig_lengths]
	feats = []

	def place(length):
		ok = [i for i, s in enumerate(seqs) if len(s) >= length + 2 * margin]
		if not ok:
			raise ValueError("no contig long enough for a feature of %d bp" % length)
		for _ in range(1000):
			ci = ok[rng.weighted_index([len(seqs[i]) - length - 2 * margin + 1 for i in ok])]
			st = margin + rng.below(len(seqs[ci]) - length - 2 * margin + 1)
			if all(not (f.contig == names[ci] and f.start < st + length + margin and st - margin < f.end) for f in feats):
				return ci, st
		raise ValueError("cannot place a %d bp feature without overlaps; genome too small" % length)

	for fam in range(repeat_families):
		unit = rng.dna(rng.randint(repeat_len[0], repeat_len[1]))
		for k in range(repeat_copies):
			ci, st = place(len(unit))
			strand = "+" if rng.below(2) == 0 else "-"
			copy = mutate(rng, unit, repeat_divergence)
			if strand == "-":
				copy = revcomp(copy)
			if soft_mask_repeats:
				copy = copy.lower()
			seqs[ci][st:st + len(copy)] = list(copy)
			feats.append(Feature(names[ci], st, st + len(copy), "repeat", (fam, k, strand)))
	for _ in range(low_complexity):
		length = rng.randint(low_complexity_len[0], low_complexity_len[1])
		period = rng.randint(1, 3)
		unit = rng.dna(period)
		if period > 1 and len(set(unit)) == 1:
			unit = unit[0] + _other_base(rng, unit[0]) + unit[2:]
		tract = (unit * (length // period + 1))[:length]
		ci, st = place(length)
		seqs[ci][st:st + length] = list(tract)
		feats.append(Feature(names[ci], st, st + length, "lowcomp", unit))
	for _ in range(n_runs):
		length = rng.randint(n_run_len[0], n_run_len[1])
		ci, st = place(length)
		seqs[ci][st:st + length] = ["N"] * length
		feats.append(Feature(names[ci], st, st + length, "N"))
	contigs = [(names[i], "".join(s)) for i, s in enumerate(seqs)]
	for k in range(short_contigs):
		contigs.append(("short%d" % (k + 1), rng.dna(rng.randint(short_contig_len[0], short_contig_len[1]))))
	feats.sort(key=lambda f: (names.index(f.contig), f.start))
	return Genome(contigs, feats)


# ---------------------------------------------------------------------------
# Read simulation
# ---------------------------------------------------------------------------


class Segment(object):
	"""Where the read bases [qstart, qend) (read orientation) came from."""
	__slots__ = ("contig", "pos", "end", "strand", "qstart", "qend", "cigar", "edits")

	def __init__(self, contig, pos, end, strand, qstart, qend, cigar, edits):
		self.contig, self.pos, self.end, self.strand = contig, pos, end, strand
		self.qstart, self.qend, self.cigar, self.edits = qstart, qend, cigar, edits

	def __repr__(self):
		return "Segment(%s:%d-%d%s q=%d-%d %s)" % (self.contig, self.pos, self.end, self.strand, self.qstart, self.qend, self.cigar)


class SimRead(object):
	"""A simulated read and its truth.

	seq/qual are as sequenced (5'->3'). For a single-segment read, contig,
	pos (0-based leftmost reference base), end (exclusive), strand and cigar
	(reference orientation, no clipping) describe the true alignment; edits is
	a list of (ref_pos, op, ref_bases, read_bases) in reference orientation
	with op 'X' (substitution), 'I' (insertion before ref_pos) or 'D'. A
	chimeric read has several segments; contig/pos/... then describe the first.
	"""

	def __init__(self, name, seq, qual, segments, comment=None):
		self.name, self.seq, self.qual, self.segments, self.comment = name, seq, qual, segments, comment

	@property
	def contig(self):
		return self.segments[0].contig

	@property
	def pos(self):
		return self.segments[0].pos

	@property
	def end(self):
		return self.segments[0].end

	@property
	def strand(self):
		return self.segments[0].strand

	@property
	def cigar(self):
		return self.segments[0].cigar

	@property
	def edits(self):
		return self.segments[0].edits

	@property
	def is_chimeric(self):
		return len(self.segments) > 1

	@property
	def fwd_seq(self):
		"""The read in reference orientation (what SAM SEQ shows for a single-segment read)."""
		return self.seq if self.strand == "+" else revcomp(self.seq)

	@property
	def nm(self):
		"""Number of simulated edits: substitutions plus indel lengths.

		Bases copied from a reference N run are not edits, so a SAM NM (which
		counts every N as a mismatch) can exceed this for such reads.
		"""
		n = 0
		for seg in self.segments:
			for _, op, rb, qb in seg.edits:
				n += 1 if op == "X" else len(rb) + len(qb)
		return n

	def __len__(self):
		return len(self.seq)

	def __repr__(self):
		return "SimRead(%s, %d bp, %r)" % (self.name, len(self.seq), self.segments)


class SimPair(object):
	"""A read pair. For FR pairs, contig/frag_start/frag_end/strand give the fragment."""

	def __init__(self, name, r1, r2, contig=None, frag_start=None, frag_end=None, strand=None):
		self.name, self.r1, self.r2 = name, r1, r2
		self.contig, self.frag_start, self.frag_end, self.strand = contig, frag_start, frag_end, strand

	@property
	def insert(self):
		return None if self.frag_start is None else self.frag_end - self.frag_start

	def __iter__(self):
		return iter((self.r1, self.r2))

	def __repr__(self):
		return "SimPair(%s, %r, %r)" % (self.name, self.r1, self.r2)


def _cigar_str(ops):
	merged = []
	for n, op in ops:
		if n <= 0:
			continue
		if merged and merged[-1][1] == op:
			merged[-1][0] += n
		else:
			merged.append([n, op])
	return "".join("%d%s" % (n, op) for n, op in merged)


def _sim_fwd(rng, ref, start, length, sub_rate, indel_rate, max_indel, end_margin, n_rate):
	"""Walk the reference from start emitting `length` read bases with edits.

	Returns (fwd_seq, end, cigar, edits) or None if the walk ran off the contig.
	Indels stay at least end_margin bases from either read end and never
	touch another indel.
	"""
	out, ops, edits = [], [], []
	i = start
	last_indel = -10
	while len(out) < length:
		if i >= len(ref):
			return None
		k = len(out)
		if indel_rate > 0 and end_margin <= k < length - end_margin - max_indel and k - last_indel > 1 and rng.chance(indel_rate):
			n = rng.randint(1, max_indel)
			if rng.below(2) == 0:  # insertion before ref[i]
				ins = rng.dna(n)
				out.extend(ins)
				ops.append((n, "I"))
				edits.append((i, "I", "", ins))
			else:
				if i + n >= len(ref):
					return None
				ops.append((n, "D"))
				edits.append((i, "D", ref[i:i + n].upper(), ""))
				i += n
			last_indel = len(out)
			continue
		c = ref[i].upper()
		if c in "ACGT" and rng.chance(sub_rate):
			b = _other_base(rng, c)
			edits.append((i, "X", c, b))
			out.append(b)
		elif c in "ACGT" and rng.chance(n_rate):
			edits.append((i, "X", c, "N"))
			out.append("N")
		else:
			out.append(c if c in "ACGT" else "N")
		ops.append((1, "M"))
		i += 1
	return "".join(out), i, _cigar_str(ops), edits


def _quals(rng, n, qual):
	if qual is None:
		return None
	if qual == "random":
		return "".join(chr(33 + 2 + rng.below(39)) for _ in range(n))
	if len(qual) == 1:
		return qual * n
	raise ValueError("qual must be None, 'random' or one character")


def _pick_contig(rng, genome, need, contigs):
	names = [n for n in (contigs or genome.names) if genome.length(n) >= need]
	if not names:
		raise ValueError("no contig of at least %d bp" % need)
	return names[rng.weighted_index([genome.length(n) - need + 1 for n in names])]


def _draw_len(rng, length):
	if isinstance(length, (tuple, list)):
		return rng.randint(length[0], length[1])
	return int(length)


def _sim_one(rng, genome, length, sub_rate, indel_rate, max_indel, end_margin, n_rate, unique, contigs,
		strand=None, contig=None, start=None, unique_pad=0):
	"""One single-segment read: returns (fwd_seq, Segment) in read orientation."""
	slack = length // 10 + 2 * max_indel + 10
	for _ in range(10000):
		c = contig or _pick_contig(rng, genome, length + slack, contigs)
		ref = genome[c]
		st = start if start is not None else rng.below(len(ref) - length - slack + 1)
		r = _sim_fwd(rng, ref, st, length, sub_rate, indel_rate, max_indel, end_margin, n_rate)
		if r is None:
			if start is not None:
				raise ValueError("read at %s:%d runs off the contig" % (c, st))
			continue
		seq, en, cigar, edits = r
		if unique and not genome.is_unique(c, st, en, unique_pad):
			if start is not None:
				raise ValueError("read at %s:%d is not in unique sequence" % (c, st))
			continue
		s = strand if strand is not None else ("+" if rng.below(2) == 0 else "-")
		return seq, Segment(c, st, en, s, 0, length, cigar, edits)
	raise ValueError("could not place a read; genome too small or too repetitive")


def simulate_reads(genome, n, length=150, seed=1, sub_rate=0.01, indel_rate=0.0005, max_indel=3, end_margin=10,
		n_rate=0.0, unique=False, unique_pad=0, contigs=None, strand=None, qual="random", name_prefix="r"):
	"""n single-end reads with recorded truth (see SimRead).

	length is an int or a (min, max) range. sub_rate and indel_rate are per
	base; indel lengths are uniform in 1..max_indel. With unique=True, reads
	come only from sequence overlapping no genome feature (within unique_pad).
	strand None picks each read's strand at random. qual is 'random', a single
	character, or None for FASTA.
	"""
	rng = Rng(seed)
	out = []
	for i in range(n):
		L = _draw_len(rng, length)
		fwd, seg = _sim_one(rng, genome, L, sub_rate, indel_rate, max_indel, end_margin, n_rate, unique, contigs,
			strand=strand, unique_pad=unique_pad)
		seq = fwd if seg.strand == "+" else revcomp(fwd)
		out.append(SimRead("%s%d" % (name_prefix, i), seq, _quals(rng, len(seq), qual), [seg]))
	return out


def simulate_long_reads(genome, n, length=(1000, 5000), seed=1, sub_rate=0.005, indel_rate=0.002, max_indel=5, **kw):
	"""Accurate long reads (HiFi-like): simulate_reads with longer lengths and more indels."""
	kw.setdefault("name_prefix", "lr")
	return simulate_reads(genome, n, length=length, seed=seed, sub_rate=sub_rate, indel_rate=indel_rate,
		max_indel=max_indel, **kw)


def simulate_pairs(genome, n, read_length=150, insert_mean=300, insert_sd=30, seed=1, sub_rate=0.01,
		indel_rate=0.0005, max_indel=3, end_margin=10, n_rate=0.0, unique=False, unique_pad=0, contigs=None,
		qual="random", name_prefix="p", min_insert=None):
	"""n FR read pairs. Inserts are ~normal(insert_mean, insert_sd), at least min_insert (default read_length).

	For a '+' fragment, r1 is forward from frag_start and r2 is reverse
	ending near frag_end; for a '-' fragment the roles swap. frag_end is the
	end of the reverse read's true reference span, so frag_end - frag_start is
	the template length a correct aligner reports for an unclipped pair.
	Both reads carry the pair's name (no /1 /2; the writers can add them).
	"""
	rng = Rng(seed)
	lo = read_length if min_insert is None else min_insert
	out = []
	for i in range(n):
		L1 = _draw_len(rng, read_length)
		L2 = _draw_len(rng, read_length)
		for _ in range(10000):
			ins = max(lo, max(L1, L2), rng.normal_int(insert_mean, insert_sd))
			slack = ins // 10 + 2 * max_indel + 10
			c = _pick_contig(rng, genome, ins + slack, contigs)
			fs = rng.below(genome.length(c) - ins - slack + 1)
			fstrand = "+" if rng.below(2) == 0 else "-"
			ref = genome[c]
			a = _sim_fwd(rng, ref, fs, L1 if fstrand == "+" else L2, sub_rate, indel_rate, max_indel, end_margin, n_rate)
			b_start = fs + ins - (L2 if fstrand == "+" else L1)
			b = _sim_fwd(rng, ref, b_start, L2 if fstrand == "+" else L1, sub_rate, indel_rate, max_indel, end_margin, n_rate)
			if a is None or b is None:
				continue
			if unique and not (genome.is_unique(c, fs, a[1], unique_pad) and genome.is_unique(c, b_start, b[1], unique_pad)):
				continue
			left = Segment(c, fs, a[1], "+", 0, len(a[0]), a[2], a[3])
			right = Segment(c, b_start, b[1], "-", 0, len(b[0]), b[2], b[3])
			lseq, rseq = a[0], revcomp(b[0])
			name = "%s%d" % (name_prefix, i)
			if fstrand == "+":
				r1 = SimRead(name, lseq, _quals(rng, len(lseq), qual), [left])
				r2 = SimRead(name, rseq, _quals(rng, len(rseq), qual), [right])
			else:
				r1 = SimRead(name, rseq, _quals(rng, len(rseq), qual), [right])
				r2 = SimRead(name, lseq, _quals(rng, len(lseq), qual), [left])
			out.append(SimPair(name, r1, r2, c, fs, b[1], fstrand))
			break
		else:
			raise ValueError("could not place a pair; genome too small")
	return out


def simulate_hic_pairs(genome, n, read_length=150, seed=1, min_distance=20000, inter_frac=0.5, chimeric_frac=0.0,
		junction_min=40, sub_rate=0.005, indel_rate=0.0, max_indel=3, unique=True, unique_pad=0, contigs=None,
		qual="random", name_prefix="h"):
	"""n Hi-C-like pairs: each mate from an independent locus.

	A fraction inter_frac of pairs put the mates on different contigs (when
	there are at least two); the others put them on one contig at least
	min_distance apart (falling back to another contig if none is long
	enough). A fraction chimeric_frac of pairs make read 1 a ligation
	chimera: its first part (>= junction_min bp) comes from its own locus,
	the rest from a third locus; SimRead.segments records both parts in read
	orientation.
	"""
	rng = Rng(seed)
	names = [c for c in (contigs or genome.names) if genome.length(c) >= read_length * 4]
	if not names:
		raise ValueError("no contig long enough for Hi-C pairs")

	def locus(L, avoid=None):
		for _ in range(10000):
			fwd, seg = _sim_one(rng, genome, L, sub_rate, indel_rate, max_indel, 10, 0.0, unique, names,
				unique_pad=unique_pad)
			if avoid is None:
				return fwd, seg
			c0, p0, inter = avoid
			if inter and seg.contig != c0:
				return fwd, seg
			if not inter and seg.contig == c0 and abs(seg.pos - p0) >= min_distance:
				return fwd, seg
		raise ValueError("cannot place a Hi-C mate; genome too small")

	out = []
	for i in range(n):
		name = "%s%d" % (name_prefix, i)
		fwd1, s1 = locus(read_length)
		long_enough = genome.length(s1.contig) >= 2 * min_distance + read_length
		inter = len(names) > 1 and (rng.chance(inter_frac) or not long_enough)
		fwd2, s2 = locus(read_length, (s1.contig, s1.pos, inter))
		reads = []
		for fwd, seg in ((fwd1, s1), (fwd2, s2)):
			seq = fwd if seg.strand == "+" else revcomp(fwd)
			reads.append([seq, [seg]])
		if rng.chance(chimeric_frac):
			a = rng.randint(junction_min, read_length - junction_min)
			# keep the first `a` bases of read 1 from its locus; take the rest from a new locus
			seq1 = reads[0][0][:a]
			b_len = read_length - a
			fwdb, sb = locus(b_len, (s1.contig, s1.pos, True if len(names) > 1 else False))
			seqb = fwdb if sb.strand == "+" else revcomp(fwdb)
			first = _trim_segment(genome, s1, a)
			sb.qstart, sb.qend = a, read_length
			reads[0] = [seq1 + seqb, [first, sb]]
		r1 = SimRead(name, reads[0][0], _quals(rng, len(reads[0][0]), qual), reads[0][1])
		r2 = SimRead(name, reads[1][0], _quals(rng, len(reads[1][0]), qual), reads[1][1])
		out.append(SimPair(name, r1, r2))
	return out


def _trim_segment(genome, seg, keep):
	"""The segment for the first `keep` read bases of an edit-free-at-the-cut segment."""
	ops = parse_cigar(seg.cigar)
	if seg.strand == "-":
		ops = ops[::-1]
	q = t = 0
	out = []
	for n, op in ops:
		if q >= keep:
			break
		if op in "MI":
			take = min(n, keep - q)
			out.append((take, op))
			q += take
			if op == "M":
				t += take
		else:
			out.append((n, op))
			t += n
	if seg.strand == "-":
		out = out[::-1]
		pos, end = seg.end - t, seg.end
	else:
		pos, end = seg.pos, seg.pos + t
	edits = [e for e in seg.edits if pos <= e[0] < end]
	return Segment(seg.contig, pos, end, seg.strand, 0, keep, _cigar_str(out), edits)


# ---------------------------------------------------------------------------
# FASTA/FASTQ I/O
# ---------------------------------------------------------------------------


def _open_out(path, gz):
	if gz is None:
		gz = path.endswith(".gz")
	if gz:
		raw = open(path, "wb")
		return gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0), raw
	return open(path, "wb"), None


def _write_text(path, text, gz=None):
	f, raw = _open_out(path, gz)
	try:
		f.write(text.encode())
	finally:
		f.close()
		if raw is not None:
			raw.close()


def _as_record(r):
	"""(name, seq, qual, comment) from a SimRead, a tuple or a dict."""
	if isinstance(r, SimRead):
		return r.name, r.seq, r.qual, r.comment
	if isinstance(r, dict):
		return r["name"], r["seq"], r.get("qual"), r.get("comment")
	r = tuple(r)
	return (r + (None, None))[:4]


def fastx_text(records, fmt="fastq", suffix=None, line_width=None):
	"""FASTA or FASTQ text. suffix is appended to every name (e.g. '/1')."""
	out = []
	for r in records:
		name, seq, qual, comment = _as_record(r)
		if suffix:
			name += suffix
		head = name if not comment else "%s %s" % (name, comment)
		if fmt == "fasta":
			out.append(">%s\n" % head)
			if line_width:
				for i in range(0, len(seq), line_width):
					out.append(seq[i:i + line_width] + "\n")
				if not seq:
					out.append("\n")
			else:
				out.append(seq + "\n")
		elif fmt == "fastq":
			if qual is None:
				qual = "I" * len(seq)
			if len(qual) != len(seq):
				raise ValueError("quality length differs from sequence length for %s" % name)
			out.append("@%s\n%s\n+\n%s\n" % (head, seq, qual))
		else:
			raise ValueError("fmt must be 'fasta' or 'fastq'")
	return "".join(out)


def _fmt_from_path(path, fmt):
	if fmt:
		return fmt
	p = path[:-3] if path.endswith(".gz") else path
	return "fasta" if p.endswith((".fa", ".fasta", ".fna")) else "fastq"


def write_fastx(path, records, fmt=None, gz=None, suffix=None, line_width=None):
	"""Write records; the format and compression follow the extension unless given."""
	_write_text(path, fastx_text(records, _fmt_from_path(path, fmt), suffix, line_width), gz)
	return path


def write_fasta(path, records, gz=None, line_width=None):
	return write_fastx(path, records, "fasta", gz, line_width=line_width)


def write_fastq(path, records, gz=None):
	return write_fastx(path, records, "fastq", gz)


def write_pairs(path1, path2, pairs, fmt=None, gz=None, suffix=False):
	"""Read 1 and read 2 into two files; suffix=True appends /1 and /2 to names."""
	write_fastx(path1, [p.r1 for p in pairs], fmt, gz, "/1" if suffix else None)
	write_fastx(path2, [p.r2 for p in pairs], fmt, gz, "/2" if suffix else None)
	return path1, path2


def write_interleaved(path, pairs, fmt=None, gz=None, suffix=False):
	"""Both mates of each pair, adjacent, in one file."""
	fmt = _fmt_from_path(path, fmt)
	parts = []
	for p in pairs:
		parts.append(fastx_text([p.r1], fmt, "/1" if suffix else None))
		parts.append(fastx_text([p.r2], fmt, "/2" if suffix else None))
	_write_text(path, "".join(parts), gz)
	return path


def _open_in(path):
	with open(path, "rb") as f:
		magic = f.read(2)
	if magic == b"\x1f\x8b":
		return io.TextIOWrapper(gzip.open(path, "rb"), encoding="ascii", errors="surrogateescape")
	return open(path, "r", encoding="ascii", errors="surrogateescape")


def read_fastx(path):
	"""[(name, seq, qual, comment)] from FASTA or FASTQ, plain or gzip; qual None for FASTA."""
	out = []
	with _open_in(path) as f:
		lines = f.read().split("\n")
	i = 0
	while i < len(lines):
		line = lines[i]
		if not line:
			i += 1
			continue
		if line[0] == ">":
			head = line[1:].split(None, 1)
			seq = []
			i += 1
			while i < len(lines) and not lines[i].startswith(">"):
				seq.append(lines[i].strip())
				i += 1
			out.append((head[0] if head else "", "".join(seq), None, head[1] if len(head) > 1 else None))
		elif line[0] == "@":
			head = line[1:].split(None, 1)
			seq = lines[i + 1]
			qual = lines[i + 3]
			out.append((head[0] if head else "", seq, qual, head[1] if len(head) > 1 else None))
			i += 4
		else:
			raise ValueError("%s: unexpected line %r" % (path, line[:50]))
	return out


def read_fasta(path):
	"""{name: sequence} from a FASTA file, plain or gzip, in file order (dict keeps order)."""
	return dict((r[0], r[1]) for r in read_fastx(path))


# ---------------------------------------------------------------------------
# Running minibwa
# ---------------------------------------------------------------------------


def _no_core():
	try:
		import resource
	except ImportError:
		return
	if hasattr(resource, "RLIMIT_CORE"):
		try:
			resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
		except (ValueError, OSError):
			pass


class Result(object):
	"""Outcome of one process: argv, returncode, stdout/stderr (str) and the raw bytes."""

	def __init__(self, argv, returncode, stdout_bytes, stderr_bytes):
		self.argv, self.returncode = argv, returncode
		self.stdout_bytes, self.stderr_bytes = stdout_bytes, stderr_bytes
		self.stdout = stdout_bytes.decode("utf-8", "surrogateescape")
		self.stderr = stderr_bytes.decode("utf-8", "surrogateescape")

	def check(self):
		if self.returncode != 0:
			raise AssertionError("%s exited with %d\nstderr:\n%s" % (" ".join(self.argv), self.returncode, self.stderr[-4000:]))
		return self

	def sam(self, strict=True):
		return parse_sam(self.stdout, strict)

	def paf(self):
		return parse_paf(self.stdout)

	def __repr__(self):
		return "Result(%r, returncode=%d)" % (" ".join(self.argv), self.returncode)


def run_minibwa(args, binary=None, timeout=None, stdin=None, env=None, cwd=None, check=False):
	"""Run `minibwa <args>`; stdin may be bytes, str or None. Core dumps are disabled."""
	argv = [binary or minibwa_build.minibwa_bin()] + [str(a) for a in args]
	if isinstance(stdin, str):
		stdin = stdin.encode()
	full_env = dict(os.environ)
	full_env.setdefault("LC_ALL", "C")
	if env:
		full_env.update(env)
	kw = {}
	if os.name == "posix":
		kw["preexec_fn"] = _no_core
	try:
		p = subprocess.run(argv, input=stdin, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
			timeout=timeout or run_timeout(), env=full_env, cwd=cwd, **kw)
	except subprocess.TimeoutExpired as e:
		raise AssertionError("%s timed out after %s s" % (" ".join(argv), e.timeout))
	r = Result(argv, p.returncode, p.stdout, p.stderr)
	return r.check() if check else r


def _has_thread_opt(args):
	for a in args:
		a = str(a)
		if a == "-t" or (a.startswith("-t") and a[2:].isdigit()):
			return True
	return False


def run_map(index, reads, args=(), threads=None, binary=None, timeout=None, stdin=None, env=None, check=True, subcommand="map"):
	"""Run `minibwa map [args] -t N <index> <reads...>`.

	reads is a path or a list of paths ('-' reads stdin). The thread count is
	`threads` if given, else whatever `args` sets with a separate -t/-tN, else
	$MINIBWA_TEST_THREADS. check=True raises on a non-zero exit.
	"""
	args = [str(a) for a in args]
	if isinstance(reads, (str, bytes)) or hasattr(reads, "__fspath__"):
		reads = [reads]
	if threads is not None:
		if _has_thread_opt(args):
			raise ValueError("give the thread count either in args or as threads=, not both")
		args = args + ["-t", str(threads)]
	elif not _has_thread_opt(args):
		args = args + ["-t", str(test_threads())]
	argv = [subcommand] + args + [str(index)] + [str(r) for r in reads]
	return run_minibwa(argv, binary=binary, timeout=timeout, stdin=stdin, env=env, check=check)


_INDEX_CACHE = {}


def index_genome(fa, tmpdir, meth=False, threads=1, binary=None, extra_args=()):
	"""Index a FASTA (path or Genome) with `minibwa index`; return the index prefix.

	Indexes are cached per (genome bytes, --meth, binary, extra args) for the
	life of the Python process, so tmpdir should outlive the test session
	(use the session-scoped `index_dir` fixture). threads sets `index -t`
	and OMP_NUM_THREADS.
	"""
	binary = binary or minibwa_build.minibwa_bin()
	tmpdir = str(tmpdir)
	if isinstance(fa, Genome):
		data = fa.fasta_text().encode()
		gh = hashlib.sha256(data).hexdigest()[:16]
		path = os.path.join(tmpdir, "genome-%s.fa" % gh)
		if not os.path.exists(path):
			with open(path, "wb") as f:
				f.write(data)
		fa = path
	else:
		fa = str(fa)
		with open(fa, "rb") as f:
			gh = hashlib.sha256(f.read()).hexdigest()[:16]
	key = (gh, bool(meth), binary, tuple(str(a) for a in extra_args))
	if key in _INDEX_CACHE and os.path.exists(_INDEX_CACHE[key] + ".l2b"):
		return _INDEX_CACHE[key]
	kh = hashlib.sha256(repr(key).encode()).hexdigest()[:12]
	d = os.path.join(tmpdir, "idx-%s" % kh)
	os.makedirs(d, exist_ok=True)
	prefix = os.path.join(d, "ref")
	argv = ["index", "-t", str(threads)] + (["--meth"] if meth else []) + [str(a) for a in extra_args] + [fa, prefix]
	run_minibwa(argv, binary=binary, env={"OMP_NUM_THREADS": str(threads)}, check=True)
	for ext in (".l2b", ".mbw") + ((".meth.mbw",) if meth else ()):
		if not os.path.exists(prefix + ext):
			raise AssertionError("minibwa index did not write %s" % (prefix + ext))
	_INDEX_CACHE[key] = prefix
	return prefix


# ---------------------------------------------------------------------------
# SAM parsing
# ---------------------------------------------------------------------------


class SamFormatError(ValueError):
	pass


_TAG_RE = re.compile(r"^([A-Za-z][A-Za-z0-9]):([AifZHB]):(.*)$", re.S)
_INT_RE = re.compile(r"^[-+]?[0-9]+$")
_FLOAT_RE = re.compile(r"^[-+]?([0-9]*\.?[0-9]+|[0-9]+\.)([eE][-+]?[0-9]+)?$")
_CIGAR_RE = re.compile(r"^(\*|([0-9]+[MIDNSHPX=])+)$")
_B_SUB = {"c": (-128, 127), "C": (0, 255), "s": (-32768, 32767), "S": (0, 65535),
	"i": (-2 ** 31, 2 ** 31 - 1), "I": (0, 2 ** 32 - 1), "f": None}


def parse_tag(text):
	"""'XX:T:value' -> (name, type, typed value). B arrays give (subtype, [values])."""
	m = _TAG_RE.match(text)
	if not m:
		raise SamFormatError("malformed tag %r" % text)
	name, typ, val = m.groups()
	if typ == "i":
		if not _INT_RE.match(val):
			raise SamFormatError("bad integer in tag %r" % text)
		v = int(val)
	elif typ == "f":
		if not _FLOAT_RE.match(val) and val not in ("inf", "-inf", "nan", "NaN"):
			raise SamFormatError("bad float in tag %r" % text)
		v = float(val)
	elif typ == "A":
		if len(val) != 1 or not ("!" <= val <= "~"):
			raise SamFormatError("bad character in tag %r" % text)
		v = val
	elif typ == "Z":
		if any(not (" " <= c <= "~") for c in val):
			raise SamFormatError("non-printable character in tag %r" % text)
		v = val
	elif typ == "H":
		if len(val) % 2 or not re.match(r"^[0-9A-F]*$", val):
			raise SamFormatError("bad hex in tag %r" % text)
		v = val
	else:
		parts = val.split(",")
		sub = parts[0]
		if sub not in _B_SUB:
			raise SamFormatError("bad B subtype in tag %r" % text)
		vals = []
		for x in parts[1:]:
			if sub == "f":
				vals.append(float(x))
			else:
				if not _INT_RE.match(x):
					raise SamFormatError("bad B array value in tag %r" % text)
				iv = int(x)
				lo, hi = _B_SUB[sub]
				if not (lo <= iv <= hi):
					raise SamFormatError("B array value out of range in tag %r" % text)
				vals.append(iv)
		v = (sub, vals)
	return name, typ, v


class SamRecord(object):
	"""One alignment line. pos is 1-based as in the file; pos0 is 0-based."""
	__slots__ = ("qname", "flag", "rname", "pos", "mapq", "cigar", "rnext", "pnext", "tlen", "seq", "qual",
		"tags", "tag_types", "tag_order", "line")

	def __init__(self, line, strict=True):
		self.line = line
		f = line.split("\t")
		if len(f) < 11:
			raise SamFormatError("fewer than 11 fields: %r" % line[:200])
		self.qname, flag, self.rname, pos, mapq, self.cigar, self.rnext, pnext, tlen, self.seq, self.qual = f[:11]
		try:
			self.flag, self.pos, self.mapq, self.pnext, self.tlen = int(flag), int(pos), int(mapq), int(pnext), int(tlen)
		except ValueError:
			raise SamFormatError("non-integer FLAG/POS/MAPQ/PNEXT/TLEN: %r" % line[:200])
		self.tags, self.tag_types, self.tag_order = {}, {}, []
		for t in f[11:]:
			name, typ, v = parse_tag(t)
			if name in self.tags:
				raise SamFormatError("duplicate tag %s in %r" % (name, line[:200]))
			self.tags[name], self.tag_types[name] = v, typ
			self.tag_order.append(name)
		if strict:
			self._validate()

	def _validate(self):
		if not re.match(r"^[!-?A-~]{1,254}$", self.qname) and self.qname != "*":
			raise SamFormatError("bad QNAME %r" % self.qname)
		if not (0 <= self.flag < 1 << 16):
			raise SamFormatError("FLAG out of range in %r" % self.line[:200])
		if not (0 <= self.pos < 1 << 31) or not (0 <= self.pnext < 1 << 31):
			raise SamFormatError("POS/PNEXT out of range in %r" % self.line[:200])
		if not (0 <= self.mapq <= 255):
			raise SamFormatError("MAPQ out of range in %r" % self.line[:200])
		if not _CIGAR_RE.match(self.cigar):
			raise SamFormatError("bad CIGAR %r" % self.cigar)
		if self.seq != "*" and not re.match(r"^[A-Za-z=.]+$", self.seq):
			raise SamFormatError("bad SEQ in %r" % self.line[:200])
		if self.qual != "*":
			if not re.match(r"^[!-~]+$", self.qual):
				raise SamFormatError("bad QUAL in %r" % self.line[:200])
			if self.seq != "*" and len(self.qual) != len(self.seq):
				raise SamFormatError("QUAL and SEQ lengths differ in %r" % self.line[:200])
		if self.cigar != "*" and self.seq != "*" and cigar_query_span(self.cigar) != len(self.seq):
			raise SamFormatError("CIGAR query length %d != SEQ length %d in %r" % (cigar_query_span(self.cigar), len(self.seq), self.line[:200]))

	def tag(self, name, default=None):
		return self.tags.get(name, default)

	@property
	def is_paired(self):
		return bool(self.flag & 0x1)

	@property
	def is_proper_pair(self):
		return bool(self.flag & 0x2)

	@property
	def is_unmapped(self):
		return bool(self.flag & 0x4)

	@property
	def mate_unmapped(self):
		return bool(self.flag & 0x8)

	@property
	def is_reverse(self):
		return bool(self.flag & 0x10)

	@property
	def mate_reverse(self):
		return bool(self.flag & 0x20)

	@property
	def is_read1(self):
		return bool(self.flag & 0x40)

	@property
	def is_read2(self):
		return bool(self.flag & 0x80)

	@property
	def is_secondary(self):
		return bool(self.flag & 0x100)

	@property
	def is_supplementary(self):
		return bool(self.flag & 0x800)

	@property
	def is_primary(self):
		return not (self.flag & 0x900)

	@property
	def segment(self):
		"""1 for read 1, 2 for read 2, 0 for an unpaired read."""
		return 1 if self.flag & 0x40 else 2 if self.flag & 0x80 else 0

	@property
	def strand(self):
		return "-" if self.flag & 0x10 else "+"

	@property
	def pos0(self):
		return self.pos - 1

	@property
	def cigar_ops(self):
		return parse_cigar(self.cigar)

	@property
	def ref_span(self):
		return cigar_ref_span(self.cigar)

	@property
	def ref_end0(self):
		"""0-based exclusive end on the reference."""
		return self.pos0 + self.ref_span

	@property
	def left_clip(self):
		"""Soft plus hard clipping at the leftmost (reference-orientation) end."""
		n = 0
		for k, op in self.cigar_ops:
			if op not in "SH":
				break
			n += k
		return n

	@property
	def right_clip(self):
		n = 0
		for k, op in reversed(self.cigar_ops):
			if op not in "SH":
				break
			n += k
		return n

	@property
	def unclipped_start0(self):
		return self.pos0 - self.left_clip

	@property
	def unclipped_end0(self):
		return self.ref_end0 + self.right_clip

	@property
	def query_length(self):
		"""Read length implied by the CIGAR, hard clips included."""
		return cigar_query_span(self.cigar, include_hard=True)

	def __repr__(self):
		return "SamRecord(%s)" % self.line[:200]


class SamFile(object):
	"""Header lines and records. hd, sq, rg and pg are parsed views of the header."""

	def __init__(self, header, records):
		self.header, self.records = header, records
		self.hd, self.sq, self.rg, self.pg, self.co = {}, [], [], [], []
		for h in header:
			code = h[1:3]
			if code == "CO":
				self.co.append(h[4:])
				continue
			fields = {}
			for x in h.split("\t")[1:]:
				if len(x) < 3 or x[2] != ":":
					raise SamFormatError("malformed header field %r in %r" % (x, h))
				if x[:2] in fields:
					raise SamFormatError("duplicate header field %s in %r" % (x[:2], h))
				fields[x[:2]] = x[3:]
			if code == "HD":
				self.hd = fields
			elif code == "SQ":
				self.sq.append((fields["SN"], int(fields["LN"])))
			elif code == "RG":
				self.rg.append(fields)
			elif code == "PG":
				self.pg.append(fields)

	def __iter__(self):
		return iter(self.records)

	def __len__(self):
		return len(self.records)

	@property
	def contig_lengths(self):
		return dict(self.sq)

	def by_qname(self):
		"""{qname: [records in file order]}."""
		d = {}
		for r in self.records:
			d.setdefault(r.qname, []).append(r)
		return d

	def primaries(self):
		"""{(qname, segment): the primary record}; raises if a segment has two primaries."""
		d = {}
		for r in self.records:
			if r.is_primary:
				k = (r.qname, r.segment)
				if k in d:
					raise SamFormatError("two primary records for %s segment %d" % k)
				d[k] = r
		return d


def parse_sam(text, strict=True):
	"""Parse SAM text. strict=True also validates field syntax per the SAM spec."""
	header, records = [], []
	for line in text.split("\n"):
		if not line:
			continue
		if line[0] == "@":
			if records:
				raise SamFormatError("header line after records: %r" % line)
			if not re.match(r"^@[A-Za-z][A-Za-z](\t|$)", line):
				raise SamFormatError("malformed header line %r" % line)
			header.append(line)
		else:
			records.append(SamRecord(line, strict))
	if text and not text.endswith("\n"):
		raise SamFormatError("output does not end with a newline")
	return SamFile(header, records)


def strip_pg(sam_text):
	"""SAM text without @PG lines (they carry the command line and the version)."""
	return "".join(l for l in sam_text.splitlines(True) if not l.startswith("@PG\t"))


def sam_body(sam_text):
	"""SAM text without any header line."""
	return "".join(l for l in sam_text.splitlines(True) if not l.startswith("@"))


# ---------------------------------------------------------------------------
# PAF parsing
# ---------------------------------------------------------------------------


class PafRecord(object):
	"""One PAF line. For an unmapped query the alignment fields are None."""
	__slots__ = ("qname", "qlen", "qstart", "qend", "strand", "tname", "tlen", "tstart", "tend", "nmatch", "alen",
		"mapq", "tags", "tag_types", "line")

	def __init__(self, line):
		self.line = line
		f = line.split("\t")
		if len(f) < 12:
			raise SamFormatError("PAF line with fewer than 12 fields: %r" % line[:200])
		self.qname = f[0]

		def num(x):
			if x == "*":
				return None
			if not _INT_RE.match(x):
				raise SamFormatError("non-integer PAF field %r in %r" % (x, line[:200]))
			return int(x)
		self.qlen, self.qstart, self.qend = num(f[1]), num(f[2]), num(f[3])
		self.strand = None if f[4] == "*" else f[4]
		self.tname = None if f[5] == "*" else f[5]
		self.tlen, self.tstart, self.tend = num(f[6]), num(f[7]), num(f[8])
		self.nmatch, self.alen, self.mapq = num(f[9]), num(f[10]), num(f[11])
		self.tags, self.tag_types = {}, {}
		for t in f[12:]:
			name, typ, v = parse_tag(t)
			if name in self.tags:
				raise SamFormatError("duplicate tag %s in %r" % (name, line[:200]))
			self.tags[name], self.tag_types[name] = v, typ

	@property
	def is_unmapped(self):
		return self.tname is None

	@property
	def is_primary(self):
		return self.tags.get("tp") == "P"

	def tag(self, name, default=None):
		return self.tags.get(name, default)

	def __repr__(self):
		return "PafRecord(%s)" % self.line[:200]


def parse_paf(text):
	out = []
	for line in text.split("\n"):
		if line:
			out.append(PafRecord(line))
	return out


# ---------------------------------------------------------------------------
# CIGAR and alignment-tag recomputation
# ---------------------------------------------------------------------------

_CIGAR_OP_RE = re.compile(r"([0-9]+)([MIDNSHPX=])")


def parse_cigar(cigar):
	"""'10M2I' -> [(10, 'M'), (2, 'I')]; '*' -> []."""
	if cigar in ("*", ""):
		return []
	ops = [(int(n), op) for n, op in _CIGAR_OP_RE.findall(cigar)]
	if "".join("%d%s" % x for x in ops) != cigar:
		raise SamFormatError("bad CIGAR %r" % cigar)
	return ops


def cigar_ref_span(cigar):
	"""Reference bases consumed: M, D, N, = and X."""
	return sum(n for n, op in parse_cigar(cigar) if op in "MDN=X")


def cigar_query_span(cigar, include_hard=False):
	"""Read bases in SEQ: M, I, S, = and X (plus H with include_hard)."""
	ops = "MIS=XH" if include_hard else "MIS=X"
	return sum(n for n, op in parse_cigar(cigar) if op in ops)


def _walk(ref, pos0, seq, cigar):
	"""Yield ('M', ref_base, read_base), ('I', read_bases), ('D', ref_bases), ('N', ref_bases)."""
	r, q = pos0, 0
	for n, op in parse_cigar(cigar):
		if op in "M=X":
			if r + n > len(ref) or q + n > len(seq):
				raise ValueError("CIGAR %s runs past the reference or the read" % cigar)
			for k in range(n):
				yield "M", ref[r + k], seq[q + k]
			r += n
			q += n
		elif op == "I":
			yield "I", seq[q:q + n]
			q += n
		elif op == "S":
			q += n
		elif op in "DN":
			if r + n > len(ref):
				raise ValueError("CIGAR %s runs past the reference" % cigar)
			yield op, ref[r:r + n]
			r += n


def _same(a, b, ambiguous_equal):
	a, b = a.upper(), b.upper()
	if a in "ACGT" and b in "ACGT":
		return a == b
	return ambiguous_equal and nt4(a) == nt4(b)


def _meth_ok(rb, qb, meth):
	rb, qb = rb.upper(), qb.upper()
	return (meth == "C2T" and rb == "C" and qb == "T") or (meth == "G2A" and rb == "G" and qb == "A")


def compute_nm(ref, pos0, seq, cigar, meth=None):
	"""Edit distance per the SAM spec: mismatches (any non-ACGT base is a mismatch) plus I and D lengths.

	meth 'C2T' or 'G2A' (reference orientation) does not count that bisulfite
	conversion as a mismatch.
	"""
	nm = 0
	for x in _walk(ref, pos0, seq, cigar):
		if x[0] == "M":
			if not _same(x[1], x[2], False) and not (meth and _meth_ok(x[1], x[2], meth)):
				nm += 1
		elif x[0] in "ID":
			nm += len(x[1])
	return nm


def compute_md(ref, pos0, seq, cigar, ambiguous_equal=False):
	"""MD:Z value in the SAM-spec form [0-9]+(([A-Z]|\\^[A-Z]+)[0-9]+)*.

	Reference bases are written uppercase. With ambiguous_equal=False (the
	htslib calmd convention) N never matches; minibwa itself treats two
	non-ACGT bases as equal.
	"""
	out, run = [], 0
	for x in _walk(ref, pos0, seq, cigar):
		if x[0] == "M":
			if _same(x[1], x[2], ambiguous_equal):
				run += 1
			else:
				out.append("%d%s" % (run, x[1].upper()))
				run = 0
		elif x[0] == "D":
			out.append("%d^%s" % (run, x[1].upper()))
			run = 0
	out.append("%d" % run)
	return "".join(out)


def compute_cs(ref, pos0, seq, cigar, long_form=False):
	"""cs:Z value (minimap2 convention, without the 'cs:Z:' prefix).

	Short form: ':N' identical run, '*xy' substitution (reference x, read y),
	'+seq' insertion, '-seq' deletion, '~xxNyy' intron; all bases lowercase.
	Bases compare by minibwa's 0-4 code, so N matches N. long_form writes
	'=SEQ' (uppercase) instead of ':N'.
	"""
	out, run = [], []

	def flush():
		if run:
			out.append("=" + "".join(run).upper() if long_form else ":%d" % len(run))
			del run[:]
	for x in _walk(ref, pos0, seq, cigar):
		if x[0] == "M":
			if nt4(x[1]) == nt4(x[2]):
				run.append(x[2])
			else:
				flush()
				out.append("*%s%s" % (x[1].lower(), x[2].lower()))
		elif x[0] == "I":
			flush()
			out.append("+" + x[1].lower())
		elif x[0] == "D":
			flush()
			out.append("-" + x[1].lower())
		else:
			flush()
			s = x[1].lower()
			out.append("~%s%d%s" % (s[:2], len(s), s[-2:]))
	flush()
	return "".join(out)


def md_from_minibwa(md):
	"""The MD value minibwa writes for a spec-form MD (no trailing 0 after a final mismatch or deletion).

	The SAMtags regex [0-9]+(([A-Z]|\\^[A-Z]+)[0-9]+)* requires the final 0;
	base minibwa (cs.c, mb_write_MD) omits it. Use only to compare against
	minibwa's MD, never to excuse a different discrepancy.
	"""
	return md[:-1] if md.endswith("0") and len(md) > 1 and not md[-2].isdigit() else md


# ---------------------------------------------------------------------------
# Generic SAM invariants
# ---------------------------------------------------------------------------


def sam_invariant_errors(sam, genome, reads=None, check_nm=True):
	"""List of violated invariants in a parsed SAM file; empty if all hold.

	genome is a Genome or a {name: sequence} dict; reads (optional) is
	{qname: [read sequences in segment order]} as written to the input. Checks:
	@SQ matches the genome; every mapped record lies within its contig;
	CIGAR query length (with hard clips) equals the read length; SEQ/QUAL of
	primary records equal the input read (reverse-complemented on the reverse
	strand); one primary per read segment; NM equals the value recomputed
	from the genome. Secondary records with SEQ '*' are checked against the
	primary's sequence.
	"""
	seqs = genome.seqs if isinstance(genome, Genome) else genome
	errs = []
	if sam.sq != [(n, len(s)) for n, s in seqs.items()]:
		errs.append("@SQ lines %r differ from the genome" % (sam.sq[:5],))
	try:
		prim = sam.primaries()
	except SamFormatError as e:
		errs.append(str(e))
		prim = {}
	seen = set()
	for r in sam.records:
		seen.add((r.qname, r.segment))
	for k in seen:
		if k not in prim:
			errs.append("no primary record for %s segment %d" % k)
	fwd = {}
	for k, p in prim.items():
		if p.seq != "*":
			fwd[k] = revcomp(p.seq) if p.is_reverse else p.seq
	if reads is not None:
		for name, rs in reads.items():
			for i, s in enumerate(rs):
				seg = 0 if len(rs) == 1 else i + 1
				k = (name, seg)
				if k not in prim:
					errs.append("read %s segment %d missing from the output" % k)
				elif k in fwd and fwd[k] != s:
					errs.append("SEQ of %s segment %d differs from the input read" % k)
	for r in sam.records:
		if r.is_unmapped:
			continue
		if r.rname not in seqs:
			errs.append("%s: unknown RNAME %s" % (r.qname, r.rname))
			continue
		clen = len(seqs[r.rname])
		if r.pos < 1 or r.ref_end0 > clen:
			errs.append("%s: alignment %s:%d-%d outside contig of length %d" % (r.qname, r.rname, r.pos, r.ref_end0, clen))
			continue
		read = fwd.get((r.qname, r.segment))
		if read is None:
			continue
		if r.query_length != len(read):
			errs.append("%s: CIGAR %s implies %d bases, read has %d" % (r.qname, r.cigar, r.query_length, len(read)))
			continue
		s = read if not r.is_reverse else revcomp(read)
		ops = r.cigar_ops
		h5 = ops[0][0] if ops and ops[0][1] == "H" else 0
		h3 = ops[-1][0] if len(ops) > 1 and ops[-1][1] == "H" else 0
		s = s[h5:len(s) - h3]
		if r.seq != "*" and r.seq != s:
			errs.append("%s: SEQ is not the read in alignment orientation" % r.qname)
			continue
		if check_nm and "NM" in r.tags:
			nm = compute_nm(seqs[r.rname], r.pos0, s, r.cigar)
			if nm != r.tags["NM"]:
				errs.append("%s: NM:i:%d but the alignment has %d differences" % (r.qname, r.tags["NM"], nm))
	return errs


def pair_invariant_errors(sam):
	"""List of violated mate-field invariants for paired records; empty if all hold.

	For every record of a two-segment read, against the mate's primary record M:
	flags 0x1, 0x40/0x80 by segment; 0x8 iff M is unmapped; 0x20 iff M is
	reverse; RNEXT/PNEXT give M's RNAME/POS ('=' when equal to this RNAME);
	'*'/0 when neither has a position; MC:Z equals M's CIGAR and MQ:i M's
	MAPQ when present. The two primaries carry opposite TLENs, zero unless
	both are mapped to one contig, and agree on flag 0x2.
	"""
	errs = []
	prim = {}
	for r in sam.records:
		if r.is_primary:
			prim[(r.qname, r.segment)] = r
	for r in sam.records:
		if r.segment == 0:
			if r.flag & 0xC9:
				errs.append("%s: unpaired record with pair flags (%d)" % (r.qname, r.flag))
			continue
		if not r.is_paired:
			errs.append("%s: segment flag without 0x1 (%d)" % (r.qname, r.flag))
		if r.is_read1 and r.is_read2:
			errs.append("%s: both 0x40 and 0x80 (%d)" % (r.qname, r.flag))
		m = prim.get((r.qname, 3 - r.segment))
		if m is None:
			errs.append("%s: no primary record for the mate of segment %d" % (r.qname, r.segment))
			continue
		if r.mate_unmapped != m.is_unmapped:
			errs.append("%s/%d: flag 0x8 is %d but the mate is %smapped" % (r.qname, r.segment, r.mate_unmapped, "un" if m.is_unmapped else ""))
		if not m.is_unmapped and r.mate_reverse != m.is_reverse:
			errs.append("%s/%d: flag 0x20 disagrees with the mate's strand" % (r.qname, r.segment))
		if m.rname == "*":
			if r.rnext != "*" or r.pnext != 0:
				errs.append("%s/%d: RNEXT/PNEXT %s/%d for a mate without position" % (r.qname, r.segment, r.rnext, r.pnext))
		else:
			want = "=" if m.rname == r.rname else m.rname
			if r.rnext != want or r.pnext != m.pos:
				errs.append("%s/%d: RNEXT/PNEXT %s/%d, mate is at %s/%d" % (r.qname, r.segment, r.rnext, r.pnext, m.rname, m.pos))
		if "MC" in r.tags and r.tags["MC"] != m.cigar:
			errs.append("%s/%d: MC:Z:%s but the mate's CIGAR is %s" % (r.qname, r.segment, r.tags["MC"], m.cigar))
		if "MQ" in r.tags and r.tags["MQ"] != m.mapq:
			errs.append("%s/%d: MQ:i:%d but the mate's MAPQ is %d" % (r.qname, r.segment, r.tags["MQ"], m.mapq))
	for (q, s), r in prim.items():
		if s != 1 or (q, 2) not in prim:
			continue
		m = prim[(q, 2)]
		if r.tlen != -m.tlen:
			errs.append("%s: TLENs %d and %d are not opposite" % (q, r.tlen, m.tlen))
		if (r.is_unmapped or m.is_unmapped or r.rname != m.rname) and r.tlen != 0:
			errs.append("%s: TLEN %d without both mates on one contig" % (q, r.tlen))
		if r.is_proper_pair != m.is_proper_pair:
			errs.append("%s: only one mate has flag 0x2" % q)
	return errs


def tlen_outer(rec, mate):
	"""TLEN as leftmost to rightmost mapped base ("TLEN#2" in the SAM spec), for two mapped
	records on one contig; positive for the leftmost segment (ties: the forward one).

	The spec leaves the definition to the implementation; minibwa uses the
	5'-to-5' form (tlen_5prime). The two agree for an FR pair whose forward
	mate starts and ends no later than the reverse mate.
	"""
	left = min(rec.pos0, mate.pos0)
	right = max(rec.ref_end0, mate.ref_end0)
	n = right - left
	first = rec.pos0 < mate.pos0 or (rec.pos0 == mate.pos0 and not rec.is_reverse)
	return n if first else -n


def tlen_5prime(rec, mate):
	"""TLEN as the 5'-to-5' distance ("TLEN#1" in the SAM spec), the form minibwa writes:
	from this record's 5' mapped base to the mate's, plus one in the direction of the
	sign, 0 when the two 5' ends coincide. Both records mapped on one contig."""
	this5 = rec.ref_end0 - 1 if rec.is_reverse else rec.pos0
	mate5 = mate.ref_end0 - 1 if mate.is_reverse else mate.pos0
	d = mate5 - this5
	return d + 1 if d > 0 else d - 1 if d < 0 else 0


def placement_error(rec, read):
	"""|unclipped start - true start| for a record of a single-segment read, or None on wrong contig/strand."""
	if rec.is_unmapped or rec.rname != read.contig or rec.strand != read.strand:
		return None
	return abs(rec.unclipped_start0 - read.pos)


# ---------------------------------------------------------------------------
# Pinned outputs
# ---------------------------------------------------------------------------


def check_pinned(name, text, max_diff_lines=40):
	"""Compare text with tests/pins/<name>; MINIBWA_TEST_UPDATE_PINS=1 rewrites the pin.

	A name ending in .gz is stored gzip-compressed. Pins are generated on
	base minibwa; a test using this must say "pins the current output" in
	its docstring.
	"""
	path = os.path.join(PINS_DIR, name)
	if os.environ.get("MINIBWA_TEST_UPDATE_PINS") == "1":
		os.makedirs(os.path.dirname(path), exist_ok=True)
		_write_text(path, text, name.endswith(".gz"))
		return
	if not os.path.exists(path):
		raise AssertionError("missing pin %s; create it on base minibwa with MINIBWA_TEST_UPDATE_PINS=1" % path)
	with _open_in(path) as f:
		want = f.read()
	if want != text:
		diff = list(difflib.unified_diff(want.splitlines(), text.splitlines(), "pinned/" + name, "actual", lineterm="", n=1))
		more = "" if len(diff) <= max_diff_lines else "\n... (%d more diff lines)" % (len(diff) - max_diff_lines)
		raise AssertionError("output differs from pin %s:\n%s%s" % (name, "\n".join(diff[:max_diff_lines]), more))


# ---------------------------------------------------------------------------
# C drivers
# ---------------------------------------------------------------------------


class CResult(object):
	def __init__(self, argv, returncode, stdout, stderr):
		self.argv, self.returncode, self.stdout, self.stderr = argv, returncode, stdout, stderr
		self.failures = [l for l in stdout.splitlines() if l.startswith("not ok")]
		m = re.search(r"^# mbt: \S+ checks=(\d+) failures=(\d+)$", stdout, re.M)
		self.checks = int(m.group(1)) if m else None
		self.n_failures = int(m.group(2)) if m else None

	@property
	def passed(self):
		return self.returncode == 0 and not self.failures and self.checks is not None and self.n_failures == 0

	def report(self):
		return "$ %s\nexit %d\n--- stdout ---\n%s\n--- stderr ---\n%s" % (" ".join(self.argv), self.returncode,
			self.stdout[-6000:], self.stderr[-3000:])


def c_cases(driver):
	"""Names of the cases a driver in tests/c/<driver>.c defines."""
	exe = minibwa_build.c_driver(driver)
	p = subprocess.run([exe, "list"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=60)
	if p.returncode != 0:
		raise AssertionError("%s list failed: %s" % (driver, p.stderr.decode("utf-8", "replace")))
	return p.stdout.decode().split()


def run_c_case(driver, case, args=(), timeout=None, stdin=None):
	"""Run one case of a C driver; return a CResult. Core dumps are disabled."""
	exe = minibwa_build.c_driver(driver)
	argv = [exe, "run", case] + [str(a) for a in args]
	kw = {}
	if os.name == "posix":
		kw["preexec_fn"] = _no_core
	try:
		p = subprocess.run(argv, input=stdin, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
			timeout=timeout or run_timeout(), **kw)
	except subprocess.TimeoutExpired as e:
		raise AssertionError("%s timed out after %s s" % (" ".join(argv), e.timeout))
	return CResult(argv, p.returncode, p.stdout.decode("utf-8", "replace"), p.stderr.decode("utf-8", "replace"))


def assert_c_case(driver, case, *args, **kw):
	"""Run a C driver case and raise AssertionError with its output unless every check passed."""
	r = run_c_case(driver, case, args, **kw)
	if not r.passed:
		raise AssertionError("C case %s/%s failed\n%s" % (driver, case, r.report()))
	return r
