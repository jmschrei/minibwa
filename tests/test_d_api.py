"""Area D: the libminibwa API (minibwa.h), through the C driver tests/c/drv_dapi.c.

mb_map() and mb_map_batch() must return, for every read, the hits that
`minibwa map` writes; the hits must not depend on the thread buffer, kalloc,
the input encoding, the index loader or the sub-batch size. The driver
prints each hit in the column layout of `minibwa map -f` (plus `sp:i:` for
mb_hit_t.sam_pri and `pp:i:` for proper_pair), so the comparison covers every public mb_hit_t and
mb_extra_t field that the PAF writer reports.
"""

import os

import pytest

import helpers_d as D
import minibwa_helpers as H


class _Data(object):
	def __init__(self, **kw):
		self.__dict__.update(kw)


@pytest.fixture(scope="module")
def mixed(synth, tmp_path_factory):
	d = tmp_path_factory.mktemp("dapi")
	recs, info = D.mixed_reads(synth.genome)
	fq = D.write_records(str(d / "mixed.fq"), recs)
	fa = H.write_fasta(str(d / "mixed.fa"), [(n, s) for n, s, q, c in recs])
	return _Data(recs=recs, info=info, fq=fq, fa=fa, dir=d)


def _dump(prefix, reads, mode="one", is_meth=0, opts=()):
	"""Run the dump case; return (hit lines without the trailing sp:i/pp:i fields, full hit lines, n_hit lines)."""
	r = H.run_c_case("drv_dapi", "dump", [prefix, is_meth, mode, reads] + list(opts))
	if not r.passed:
		raise AssertionError(r.report())
	hits, full, counts = [], [], []
	for l in r.stdout.splitlines():
		if l.startswith("H\t"):
			full.append(l[2:])
			hits.append(l[2:l.rindex("\tsp:i:")])
		elif l.startswith("N\t"):
			counts.append(l[2:])
	return hits, full, counts


def _paf_hits(text):
	"""Mapped PAF lines of `minibwa map -f` (unmapped lines dropped)."""
	return [l for l in text.splitlines() if l.split("\t")[5] != "*"]


def test_idx_load_and_contigs(synth, index_of, tmp_path):
	"""mb_idx_load() and mb_idx_load_mmap() return NULL for a missing prefix and for is_meth on an index
	without .meth.mbw, load a plain and a --meth index, and report contig names and lengths from the FASTA;
	mb_idx_ctg_name()/mb_idx_ctg_len() return NULL/-1 out of range; mb_idx_destroy(NULL) is a no-op."""
	meth = index_of(synth.genome, meth=True)
	args = [synth.index, str(tmp_path / "no-such-index"), meth, len(synth.genome.names)]
	for n in synth.genome.names:
		args += [n, synth.genome.length(n)]
	H.assert_c_case("drv_dapi", "idx_load", *args)


def test_opt_init_and_presets():
	"""mb_opt_init() equals the adap preset, mb_opt_preset() resets every field, and returns -1 for an
	unknown preset; the defaults match the `minibwa map` usage text (-t 1, --outn 0, --outs 0.8, --xa 5,
	-K 100m,1g) and seed 11."""
	H.assert_c_case("drv_dapi", "opt")


@pytest.mark.parametrize("flags,cli", [((), ()), (("cs",), ("-b", "cs")), (("MD",), ("-b", "MD")), (("ds",), ("-b", "ds")),
	(("eqx",), ("--eqx",))], ids=["plain", "cs", "MD", "ds", "eqx"])
def test_mb_map_equals_minibwa_map(synth, mixed, flags, cli):
	"""Every hit mb_map() returns (one read at a time, NULL thread buffer) equals the corresponding line of
	`minibwa map -f --outn 1000000 --outs 0`, which writes every hit: coordinates, mlen/blen, MAPQ, tp, s1,
	cm, s2, NM = blen-mlen+n_ambi, AS = dp_score, ms = dp_max0, md = dp_max-dp_max2, the CIGAR decoded with
	MB_CIGAR_STR, and the cs/MD/ds string stored after cigar[n_cigar] when mb_extra_t.cs is set."""
	hits, _, counts = _dump(synth.index, mixed.fa, opts=flags)
	paf = H.run_map(synth.index, mixed.fa, ["-f", "--outn", "1000000", "--outs", "0"] + list(cli)).stdout
	assert hits == _paf_hits(paf)
	assert [c.split("\t")[0] for c in counts] == [r[0] for r in mixed.recs]
	n_unmapped = sum(1 for l in paf.splitlines() if l.split("\t")[5] == "*")
	assert n_unmapped == sum(1 for c in counts if c.endswith("\t0")) >= 3
	if flags == ("cs",):
		assert all("\tcs:Z:" in h for h in hits)


def test_api_modes_give_identical_hits(synth, mixed):
	"""mb_map() with a kalloc or malloc thread buffer, with mb_tbuf_reset() forcing a reset after every
	read, with 0/1/2/3-encoded input, and on an index loaded with mb_idx_load_mmap(), and
	mb_map_batch() (MB_F_PE cleared) with sub-batches of one read or one base, n_thread 7 and a
	thread buffer, all return the same hits (MB_F_PE, n_thread and the buffer are not used by mb_map)."""
	base, base_full, _ = _dump(synth.index, mixed.fa)
	runs = [("one_tbuf", ()), ("one_tbuf_nok", ()), ("one_reset", ()), ("one_nt4", ()), ("one", ("mmap=1",)), ("one", ("mmap=0",)),
		("one", ("-pe", "nthread=7")), ("batch", ("-pe",)), ("batch", ("-pe", "sbseq=1")), ("batch", ("-pe", "sblen=1")),
		("batch_tbuf", ("-pe", "nthread=7"))]
	for mode, opts in runs:
		_, full, _ = _dump(synth.index, mixed.fa, mode, opts=opts)
		assert full == base_full, (mode, opts)


def test_batch_pairs_adjacent_reads(synth, tmp_path):
	"""With MB_F_PE (the mb_opt_init() default) mb_map_batch() pairs reads 2i and 2i+1, the last pair
	included: mates of FR pairs (inserts within the default bounds 50-800) get proper_pair hits at their
	true positions; with MB_F_NO_PAIRING the hits equal mb_map()'s."""
	pairs = H.simulate_pairs(synth.genome, 30, seed=71, unique=True)
	fa = H.write_fasta(str(tmp_path / "il.fa"), [r for p in pairs for r in (p.r1, p.r2)])
	_, full, _ = _dump(synth.index, fa, "batch")
	byname = {}
	for l in full:
		f = l.split("\t")
		byname.setdefault(f[0], []).append(f)
	for p in pairs:
		assert 50 <= p.insert <= 800
		got = byname[p.name]
		pri = [f for f in got if f[-2] == "sp:i:1"]
		assert len(pri) == 2
		for f, r in zip(pri, (p.r1, p.r2)):
			assert f[5] == r.contig and f[4] == r.strand
			assert abs(int(f[7]) - (r.pos + (int(f[2]) if r.strand == "+" else len(r.seq) - int(f[3])))) <= 3
			assert f[-1] == "pp:i:1", p.name
	_, nopair, _ = _dump(synth.index, fa, "batch", opts=("nopair",))
	_, one, _ = _dump(synth.index, fa, "one")
	assert nopair == one


def test_eqx_cigar_against_genome(synth, mixed):
	"""With MB_F_EQX every CIGAR uses MB_CIGAR_EQ_MATCH/MB_CIGAR_X_MISMATCH ('='/'X') instead of 'M', and
	each '=' base equals the reference base and each 'X' base differs (or is not ACGT)."""
	hits, _, _ = _dump(synth.index, mixed.fa, opts=("eqx",))
	reads = dict((r[0], r[1]) for r in mixed.recs)
	n_x = 0
	for h in hits:
		f = h.split("\t")
		name, qlen, qs, qe, strand, ctg, ts = f[0], int(f[1]), int(f[2]), int(f[3]), f[4], f[5], int(f[7])
		cg = [t for t in f if t.startswith("cg:Z:")][0][5:]
		assert "M" not in cg
		q = reads[name] if strand == "+" else H.revcomp(reads[name])
		q = q[qs:qe] if strand == "+" else q[qlen - qe:qlen - qs]
		ref = synth.genome[ctg]
		i = j = 0
		for n, op in H.parse_cigar(cg):
			if op in "=X":
				for k in range(n):
					same = ref[ts + j + k].upper() == q[i + k].upper() and q[i + k].upper() in "ACGT"
					assert same == (op == "="), (name, cg, i + k)
				n_x += n if op == "X" else 0
				i += n
				j += n
			elif op == "I":
				i += n
			elif op == "D":
				j += n
		assert i == qe - qs
	assert n_x > 20


def test_batch_edge_cases(synth, mixed):
	"""mb_map_batch() returns NULL for n_seq <= 0 and for MB_F_METH on an index loaded without is_meth;
	mb_map() with mt 1 or 2 on such an index returns no hits; qname may be NULL."""
	H.assert_c_case("drv_dapi", "batch_edge", synth.index, mixed.fa)


def test_seed_breaks_ties(synth, tmp_path):
	"""mb_opt_t.seed feeds the hash that breaks ties between equally good hits: with the default seed (11)
	the SAM-primary copy of each exact-repeat read equals `minibwa map`'s; every seed picks one of the
	true copies; and different seeds do not all pick the same copies."""
	reps = [f for f in synth.genome.features if f.kind == "repeat"]
	fam0 = [f for f in reps if f.info[0] == 0]
	unit_len = min(f.end - f.start for f in fam0)
	recs = []
	for i in range(24):
		f = fam0[0]
		st = f.start + 10 + (i * 37) % (unit_len - 170)
		s = synth.genome[f.contig][st:st + 150]
		recs.append(("t%d" % i, s if f.info[2] == "+" else H.revcomp(s)))
	fa = H.write_fasta(str(tmp_path / "rep.fa"), recs)
	r = H.run_c_case("drv_dapi", "seed", [synth.index, fa, 6])
	assert r.passed, r.report()
	rows = [l.split("\t")[2:] for l in r.stdout.splitlines() if l.startswith("S\t")]
	assert len(rows) == 6 and all(len(x) == len(recs) for x in rows)
	starts = set((f.contig, f.start) for f in fam0)
	for row in rows:
		for x in row:
			ctg, ts, _ = x.split(":")
			# the copies are exact, so a hit starts at the same offset in some copy
			assert any(c == ctg and s <= int(ts) < s + unit_len for c, s in starts)
	sam = H.run_map(synth.index, fa).sam()
	cli = ["%s:%d:%s" % (p.rname, p.pos0, p.strand) for p in sorted(sam.primaries().values(), key=lambda x: int(x.qname[1:]))]
	assert rows[0] == cli
	assert len(set(tuple(x) for x in rows)) > 1


def test_mb_format_matches_minibwa_map(synth, mixed):
	"""mb_format() on mb_map() hits, with the read's name, sequence and qualities, writes the same SAM and
	PAF records as `minibwa map --outn 1000000 --outs 0` (every hit written, unmapped reads included) when
	mb_opt_t.out_s is 0 too (out_s also limits XA:Z)."""
	r = H.run_c_case("drv_dapi", "fmt", [synth.index, mixed.fq, "outs=0"])
	assert r.passed, r.report()
	api = "".join(l + "\n" for l in r.stdout.splitlines() if not l.startswith("#"))
	assert api == H.sam_body(H.run_map(synth.index, mixed.fq, ["--outn", "1000000", "--outs", "0"]).stdout)
	r = H.run_c_case("drv_dapi", "fmt", [synth.index, mixed.fq, "paf"])
	assert r.passed, r.report()
	api = "".join(l + "\n" for l in r.stdout.splitlines() if not l.startswith("#"))
	assert api == H.run_map(synth.index, mixed.fq, ["-f", "--outn", "1000000", "--outs", "0"]).stdout


def _fmt_records(prefix, fq, *opts):
	r = H.run_c_case("drv_dapi", "fmt", [prefix, fq] + list(opts))
	assert r.passed, r.report()
	return H.parse_sam("".join(l + "\n" for l in r.stdout.splitlines() if not l.startswith("#"))).records


def test_mb_format_2nd_seq(synth, mixed):
	"""MB_F_2ND_SEQ (API only) writes SEQ/QUAL on secondary records: the aligned part of the read
	(reverse-complemented and reversed on the reverse strand) with hard clips; with MB_F_SUPP_SOFT too,
	the whole read with soft clips. Other records are unchanged."""
	plain = _fmt_records(synth.index, mixed.fq)
	sec = _fmt_records(synth.index, mixed.fq, "2nd")
	soft = _fmt_records(synth.index, mixed.fq, "2nd", "soft")
	reads = dict((r[0], (r[1], r[2])) for r in mixed.recs)
	assert len(plain) == len(sec) == len(soft)
	n_sec = n_clipped = 0
	for a, b, c in zip(plain, sec, soft):
		if not a.is_secondary:
			assert a.line == b.line
			continue
		n_sec += 1
		assert a.seq == a.qual == "*"
		seq, qual = reads[a.qname]
		if a.is_reverse:
			seq, qual = H.revcomp(seq), qual[::-1]
		ops = H.parse_cigar(b.cigar)
		assert all(op != "S" for _, op in ops)
		lc, rc = b.left_clip, b.right_clip
		n_clipped += lc + rc > 0
		assert b.seq == seq[lc:len(seq) - rc] and b.qual == qual[lc:len(qual) - rc]
		assert b.cigar.replace("H", "S") == a.cigar
		assert c.cigar == a.cigar and c.seq == seq and c.qual == qual
		assert b.line.split("\t")[11:] == a.line.split("\t")[11:]
	assert n_sec >= 10 and n_clipped >= 2


@pytest.fixture(scope="module")
def bs(synth, index_of, tmp_path_factory):
	d = tmp_path_factory.mktemp("dapi_bs")
	pairs = D.bisulfite_pairs(synth.genome, 24, seed=73, sub_rate=0.005)
	se = H.write_fasta(str(d / "se.fa"), [(n, r1) for n, r1, r2, _ in pairs])
	il = H.write_fasta(str(d / "il.fa"), [x for n, r1, r2, _ in pairs for x in ((n, r1), (n, r2))])
	return _Data(index=index_of(synth.genome, meth=True), pairs=pairs, se=se, il=il)


def test_meth_api_matches_minibwa_map(bs):
	"""On a --meth index, mb_map() with mt=1 (read 1, C-to-T) equals `minibwa map --meth -f` on single-end
	reads; mb_map() with mt alternating 1/2 equals mb_map_batch() with MB_F_METH|MB_F_PE|MB_F_NO_PAIRING
	(read 2i+1 converted G-to-A) and `minibwa map --meth -P -f` on the interleaved pairs; with MB_F_PE
	cleared, mb_map_batch() treats every read as read 1."""
	hits, _, _ = _dump(bs.index, bs.se, "one", 1, ("mt=1",))
	assert len(hits) >= len(bs.pairs)
	assert hits == _paf_hits(H.run_map(bs.index, bs.se, ["--meth", "-f", "--outn", "1000000", "--outs", "0"]).stdout)
	_, par, _ = _dump(bs.index, bs.il, "one_parity", 1)
	_, bat, _ = _dump(bs.index, bs.il, "batch", 1, ("meth", "nopair"))
	assert par == bat
	cli = _paf_hits(H.run_map(bs.index, bs.il, ["--meth", "-P", "-f", "--outn", "1000000", "--outs", "0"]).stdout)
	_, _, counts = _dump(bs.index, bs.il, "one_parity", 1)
	mine, it = [], iter(par)
	for i, c in enumerate(counts):  # reads alternate r1/r2 under one name; PAF names carry /1 and /2
		for _ in range(int(c.split("\t")[1])):
			l = next(it)
			mine.append("%s/%d%s" % (c.split("\t")[0], i % 2 + 1, l[l.index("\t"):l.rindex("\tsp:i:")]))
	assert mine == cli
	_, se1, _ = _dump(bs.index, bs.il, "one", 1, ("mt=1",))
	_, se_batch, _ = _dump(bs.index, bs.il, "batch", 1, ("meth", "-pe"))
	assert se1 == se_batch
