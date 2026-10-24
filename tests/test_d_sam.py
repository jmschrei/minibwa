"""Area D: the SAM writer (format.c) and the output step of `minibwa map`.

One read set (helpers_d.mixed_reads: unique, repeat, chimeric, unmapped,
short, long, IUPAC/lowercase and N-run reads, some with comments) is mapped
once per option set and each output is checked many ways against the input
reads, the genome and the SAM specification. A second genome
(helpers_d.GradedRepeats) has repeat copies whose scores are known in
advance, for XA:Z, n2:i, --outs, --xa, --xa-ratio and --outn.
"""

import pytest

import helpers_d as D
import minibwa_helpers as H

RG = "@RG\\tID:grp1\\tSM:s1\\tPL:ILLUMINA"
CO = "@CO\\thello world"


class _Data(object):
	def __init__(self, **kw):
		self.__dict__.update(kw)


@pytest.fixture(scope="module")
def mixed(synth, tmp_path_factory):
	d = tmp_path_factory.mktemp("dsam")
	recs, info = D.mixed_reads(synth.genome)
	fq = D.write_records(str(d / "mixed.fq"), recs)
	runs = {
		"default": H.run_map(synth.index, fq),
		"opts": H.run_map(synth.index, fq, ["--outn", "3", "-b", "cs", "-y", "-R", RG, "-H", CO]),
		"Y": H.run_map(synth.index, fq, ["-Y"]),
		"u": H.run_map(synth.index, fq, ["-u"]),
	}
	return _Data(recs=recs, info=info, fq=fq, runs=runs, sams=dict((k, v.sam()) for k, v in runs.items()),
		reads=dict((r[0], r) for r in recs))


def _oriented(read, rec):
	"""(SEQ, QUAL) of the whole read in the record's orientation."""
	_, seq, qual, _ = read
	if rec.is_reverse:
		return H.revcomp(seq), qual[::-1]
	return seq, qual


def test_header(synth, mixed):
	"""@HD VN:1.6 SO:unsorted GO:query; one @SQ per contig in FASTA order with its length; the -R line with
	escapes expanded; @PG with ID, PN, VN = MB_VERSION and CL = the command line; then the -H lines."""
	res = mixed.runs["opts"]
	head = [l for l in res.stdout.split("\n") if l.startswith("@")]
	want = ["@HD\tVN:1.6\tSO:unsorted\tGO:query"]
	want += ["@SQ\tSN:%s\tLN:%d" % (n, synth.genome.length(n)) for n in synth.genome.names]
	want += ["@RG\tID:grp1\tSM:s1\tPL:ILLUMINA"]
	want += ["@PG\tID:minibwa\tPN:minibwa\tVN:%s\tCL:minibwa %s" % (D.mb_version(), " ".join(res.argv[1:]))]
	want += ["@CO\thello world"]
	assert head == want
	res0 = mixed.runs["default"]
	head0 = [l for l in res0.stdout.split("\n") if l.startswith("@")]
	assert head0 == want[:-3] + ["@PG\tID:minibwa\tPN:minibwa\tVN:%s\tCL:minibwa %s" % (D.mb_version(), " ".join(res0.argv[1:]))]


@pytest.mark.parametrize("run", ["default", "opts", "Y", "u"])
def test_sam_invariants(synth, mixed, run):
	"""Strict SAM syntax; @SQ from the genome; one primary per read; SEQ of primaries is the input read
	(reverse-complemented on the reverse strand, IUPAC codes and case kept); CIGAR (hard clips included)
	spans the read; alignments lie inside their contig; NM equals the recomputed edit distance."""
	sam = mixed.sams[run]
	reads = dict((n, [s]) for n, s, q, c in mixed.recs)
	if run == "u":
		reads = dict((n, v) for n, v in reads.items() if mixed.info[n]["cat"] != "random")
	assert H.sam_invariant_errors(sam, synth.genome, reads) == []
	assert all(r.segment == 0 and not r.flag & 0x6EB for r in sam.records)


def test_input_order(mixed):
	"""Records come out in input order (by first appearance of each read name) at any thread count."""
	want = [r[0] for r in mixed.recs]
	for k, sam in mixed.sams.items():
		seen = []
		for r in sam.records:
			if not seen or seen[-1] != r.qname:
				seen.append(r.qname)
		assert len(seen) == len(set(seen)), k
		assert seen == [n for n in want if n in set(seen)], k


def test_placement_and_strand_flag(mixed):
	"""Reads from unique sequence (plain, short, long, IUPAC, N-run) have their primary at the true contig
	and strand (flag 0x10 iff reverse) with the unclipped start within 5 bp of the truth; each part of a
	chimeric read is aligned at its true locus as the primary or a supplementary (0x800) record."""
	sam = mixed.sams["default"]
	prim = sam.primaries()
	byname = sam.by_qname()
	n = 0
	for name, d in mixed.info.items():
		if d["cat"] in ("unique", "short", "long", "iupac", "nrun", "repeat_edge"):
			rec = prim[(name, 0)]
			assert not rec.is_unmapped and rec.rname == d["contig"], rec.line[:120]
			assert rec.is_reverse == (d["strand"] == "-")
			assert abs(rec.unclipped_start0 - d["pos"]) <= 5, (rec.line[:120], d)
			n += 1
		elif d["cat"] == "chimeric":
			recs = [r for r in byname[name] if not r.is_secondary]
			assert sum(r.is_primary for r in recs) == 1
			assert all(r.is_primary or r.is_supplementary for r in recs)
			for qs, qe, ctg, pos, strand in d["parts"]:
				hits = [r for r in recs if r.rname == ctg and r.strand == strand and
					abs(r.pos0 - (pos + (r.left_clip - (qs if strand == "+" else 150 - qe)))) <= 3]
				assert len(hits) == 1, (name, d["parts"], [r.line[:80] for r in recs])
	assert n >= 50


def test_unmapped_records(mixed):
	"""An unmapped single-end read: FLAG 4, RNAME '*', POS 0, MAPQ 0, CIGAR '*', RNEXT '*', PNEXT 0,
	TLEN 0, SEQ and QUAL as in the input, no alignment tags; with -R it still carries RG:Z."""
	for run in ("default", "opts"):
		recs = [r for r in mixed.sams[run].records if mixed.info[r.qname]["cat"] == "random"]
		assert len(recs) == 3
		for r in recs:
			_, seq, qual, _ = mixed.reads[r.qname]
			f = r.line.split("\t")
			assert f[1:11] == ["4", "*", "0", "0", "*", "*", "0", "0", seq, qual]
			assert set(r.tags) <= set(["RG"])
			assert ("RG" in r.tags) == (run == "opts")


def test_u_drops_only_unmapped(mixed):
	"""-u drops exactly the unmapped records; all other bytes are unchanged."""
	body = H.sam_body(mixed.runs["default"].stdout).splitlines(True)
	kept = [l for l in body if not int(l.split("\t")[1]) & 4]
	assert len(kept) == len(body) - 3
	assert H.sam_body(mixed.runs["u"].stdout).splitlines(True) == kept


def test_seq_qual_and_clipping(mixed):
	"""SEQ/QUAL are the read reverse-complemented and reversed on the reverse strand; secondary records
	have SEQ and QUAL '*'; supplementary records are hard-clipped and carry only the aligned bases, while
	primary records soft-clip; with -Y supplementary records soft-clip and carry the whole read, and every
	other record is unchanged."""
	n_supp = n_sec = 0
	for run in ("default", "opts", "Y"):
		for r in mixed.sams[run].records:
			if r.is_unmapped:
				continue
			seq, qual = _oriented(mixed.reads[r.qname], r)
			ops = H.parse_cigar(r.cigar)
			if r.is_secondary:
				assert r.seq == "*" and r.qual == "*"
				assert all(op != "H" for _, op in ops)
				n_sec += 1
			elif r.is_supplementary and run != "Y":
				assert all(op != "S" for _, op in ops), r.line[:100]
				lc, rc = r.left_clip, r.right_clip
				assert lc + rc > 0
				assert r.seq == seq[lc:len(seq) - rc] and r.qual == qual[lc:len(qual) - rc]
				n_supp += 1
			else:
				assert all(op != "H" for _, op in ops), r.line[:100]
				assert r.seq == seq and r.qual == qual
	assert n_sec >= 8 and n_supp >= 8
	a = H.sam_body(mixed.runs["default"].stdout).splitlines()
	b = H.sam_body(mixed.runs["Y"].stdout).splitlines()
	assert len(a) == len(b)
	for x, y in zip(a, b):
		fx, fy = x.split("\t"), y.split("\t")
		if int(fx[1]) & 0x800:
			assert fy[5] == fx[5].replace("H", "S") and fx[:5] + fx[6:9] + fx[11:] == fy[:5] + fy[6:9] + fy[11:]
		else:
			assert x == y


def test_sa_tag(mixed):
	"""Every non-secondary record of a read with several non-secondary records carries SA:Z listing the
	others in output order as rname,pos,strand,CIGAR,mapQ,NM; with CIGAR = 5' clip S, min(query, ref span)
	M, the excess as I or D, 3' clip S (SAM spec). Reads with one such record carry no SA:Z. Half the
	chimeric reads carry an N in each part, which NM counts."""
	n = n_nm = 0
	for run in ("default", "Y"):
		for name, recs in mixed.sams[run].by_qname().items():
			main = [r for r in recs if not r.is_secondary and not r.is_unmapped]
			qlen = len(mixed.reads[name][1])
			for r in recs:
				if r.is_secondary or len(main) < 2:
					assert "SA" not in r.tags, r.line[:150]
					continue
				sa = D.parse_sa(r.tags["SA"])
				assert sa == D.expected_sa(main, r, qlen), r.line[:200]
				n += 1
				n_nm += sum(e[5] > 0 for e in sa)
	assert n >= 32 and n_nm >= 8


def test_tags_from_alignment(synth, mixed):
	"""For records whose CIGAR is one M operation (an ungapped alignment), AS:i is the sum of the column
	scores (match 2, mismatch -8, ambiguous base -1) and ms:i the best-scoring segment, both recomputed
	from the genome; md:i is at most ms:i and equals it for a uniquely placed read with no secondary hit."""
	n = 0
	byname = mixed.sams["opts"].by_qname()
	for r in mixed.sams["default"].records:
		if r.is_unmapped:
			continue
		assert 0 <= r.tags["md"]
		if mixed.info[r.qname]["cat"] != "long":
			assert r.tags["md"] <= r.tags["ms"]
		ops = H.parse_cigar(r.cigar)
		m = [k for k, op in ops if op == "M"]
		if len(m) != 1 or any(op in "ID" for _, op in ops) or r.is_secondary:
			continue
		seq = r.seq[r.left_clip:len(r.seq) - r.right_clip] if not r.is_supplementary else r.seq
		AS, ms = D.ungapped_scores(synth.genome[r.rname], r.pos0, seq)
		if r.left_clip + r.right_clip == 0:
			assert r.tags["AS"] == AS, r.line[:150]
		assert r.tags["ms"] == ms, r.line[:150]
		if mixed.info[r.qname]["cat"] == "unique" and r.mapq == 60 and len(byname[r.qname]) == 1:
			assert r.tags["md"] == r.tags["ms"]
		n += 1
	assert n >= 40


def test_rg_comment_and_tag_options(mixed):
	"""With -R every record (unmapped and secondary included) carries RG:Z:<ID>; with -y each record of a
	read with a FASTQ comment ends with that comment; -b cs adds cs:Z after the alignment tags. These
	options, and --outn adding secondary records, leave the other fields of every record unchanged."""
	opts = mixed.sams["opts"]
	for r in opts.records:
		assert r.tags.get("RG") == "grp1"
		comment = mixed.reads[r.qname][3]
		if comment:
			assert r.line.endswith("\t" + comment)
		else:
			assert r.tag_order[-1] in ("RG", "NM", "AS", "ms", "md", "cs", "SA", "n2", "XA")
		assert ("cs" in r.tags) == (not r.is_unmapped)
	assert all("RG" not in r.tags and "cs" not in r.tags for r in mixed.sams["default"].records)

	def strip(rec):
		f = rec.line.split("\t")
		keep = f[:11] + [t for t in f[11:] if t[:2] not in ("RG", "cs") and t != mixed.reads[rec.qname][3]]
		return "\t".join(keep)
	a = [r.line for r in mixed.sams["default"].records]
	b = [strip(r) for r in opts.records if not r.is_secondary]
	assert a == b


def test_iupac_and_lowercase_seq(mixed):
	"""SEQ keeps IUPAC codes and lower case: verbatim on the forward strand, complemented per IUPAC (R<->Y,
	K<->M, B<->V, D<->H, S, W, N unchanged) with case kept on the reverse strand; NM counts every non-ACGT
	read base."""
	n = 0
	for r in mixed.sams["default"].records:
		if mixed.info[r.qname]["cat"] != "iupac" or not r.is_primary:
			continue
		seq = mixed.reads[r.qname][1]
		assert r.seq == (H.revcomp(seq) if r.is_reverse else seq)
		assert any(c.islower() for c in r.seq) == any(c.islower() for c in seq)
		assert r.cigar == "150M" and r.tags["NM"] >= 3
		n += 1
	assert n == 4


@pytest.mark.pinned
def test_pin_mixed(mixed):
	"""Pins the current output: SAM of the mixed read set, default options and with --outn 3 -b cs -y -R -H."""
	H.check_pinned("d_mixed.sam.gz", H.strip_pg(mixed.runs["default"].stdout))
	H.check_pinned("d_mixed_opts.sam.gz", H.strip_pg(mixed.runs["opts"].stdout))


# ---------------------------------------------------------------------------
# XA:Z, n2:i, --outs, --xa, --xa-ratio and --outn on hits known in advance
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def graded(index_of, tmp_path_factory):
	g = D.GradedRepeats()
	d = tmp_path_factory.mktemp("dxa")
	fq = H.write_fastq(str(d / "g.fq"), g.reads())
	return _Data(g=g, index=index_of(g.genome), fq=fq, cache={})


def _graded_run(graded, *args):
	key = tuple(args)
	if key not in graded.cache:
		graded.cache[key] = H.run_map(graded.index, graded.fq, list(args))
	return graded.cache[key]


def _flip(strand, rec):
	return strand if not rec.is_reverse else "+-"[strand == "+"]


def _xa_want(copies, prim):
	"""XA entries of the copies for a primary record; the read of a reverse primary is reverse-complemented."""
	out = []
	for c in copies:
		s = c["strand"] if prim.qname == "q" else "+-"[c["strand"] == "+"]
		out.append((c["ctg"], s, c["pos"], "150M", c["nm"]))
	return sorted(out)


@pytest.mark.parametrize("args,outs,xa", [((), 0.8, 5), (("--outs", "0.75"), 0.75, 5), (("--xa-ratio", "0.75"), 0.75, 5),
	(("--outs", "0.625"), 0.625, 5), (("--outs", "0.625", "--xa", "4"), 0.625, 4), (("--outs", "0.625", "--xa", "3"), 0.625, 3),
	(("--xa", "0",), 0.8, 0), (("--xa", "6"), 0.8, 6)], ids=["default", "outs.75", "xa-ratio.75", "outs.625", "outs.625-xa4",
	"outs.625-xa3", "xa0", "xa6"])
def test_xa_and_n2(graded, args, outs, xa):
	"""n2:i counts the secondary hits with dp_max >= --outs times the primary's (the bound included;
	--xa-ratio is the same option); XA:Z lists exactly those hits as ctg,strandPos,CIGAR,NM when there
	are at most --xa of them; --xa 0 drops both tags. Family A's copies score 280 (primary), 260, 250,
	210 and 200; family B's read has six identical secondaries, so with --xa 5 it gets n2:i:6 and no XA."""
	g = graded.g
	sam = _graded_run(graded, *args).sam()
	prim = sam.primaries()
	assert len(sam.records) == 3
	for q in ("q", "qr"):
		r = prim[(q, 0)]
		c0 = g.copies[0]
		assert (r.rname, r.pos, r.strand) == (c0["ctg"], c0["pos"], c0["strand"] if q == "q" else "-")
		assert r.cigar == "150M" and r.tags["NM"] == c0["nm"] and r.tags["AS"] == r.tags["ms"] == c0["score"]
		assert r.tags["md"] == c0["score"] - g.copies[1]["score"]
		sec = g.secondaries_at(outs)
		if xa == 0:
			assert "n2" not in r.tags and "XA" not in r.tags
			continue
		assert r.tags["n2"] == len(sec)
		if len(sec) <= xa:
			assert sorted(D.parse_xa(r.tags["XA"])) == _xa_want(sec, r)
		else:
			assert "XA" not in r.tags
	b = prim[("b", 0)]
	assert b.tags["NM"] == 0 and b.tags["md"] == 0
	if xa == 0:
		assert "n2" not in b.tags
	else:
		assert b.tags["n2"] == 6
		if xa >= 6:
			hits = [(c, s, p) for c, s, p, cig, nm in D.parse_xa(b.tags["XA"])] + [(b.rname, b.strand, b.pos)]
			assert sorted(hits) == sorted((c["ctg"], c["strand"], c["pos"]) for c in g.copies_b)
			assert all(cig == "150M" and nm == 0 for c, s, p, cig, nm in D.parse_xa(b.tags["XA"]))
		else:
			assert "XA" not in b.tags


def test_xa_ratio_is_outs(graded):
	"""--xa-ratio and --outs set the same value: identical output apart from @PG."""
	a = _graded_run(graded, "--outs", "0.75", "--outn", "10")
	b = _graded_run(graded, "--xa-ratio", "0.75", "--outn", "10")
	assert H.strip_pg(a.stdout) == H.strip_pg(b.stdout)


@pytest.mark.parametrize("outn,outs", [(10, None), (10, "0.75"), (10, "0"), (1, None), (2, "0"), (0, "0")])
def test_outn_secondary_records(graded, outn, outs):
	"""--outn N writes up to N secondary records (FLAG 0x100, SEQ/QUAL '*', MAPQ 0) per read, the best
	first, among the hits passing --outs; each has the copy's position, strand, CIGAR, NM, AS and ms; the
	primary record is unchanged."""
	args = ["--outn", str(outn)] + (["--outs", outs] if outs else [])
	sam = _graded_run(graded, *args).sam()
	base = _graded_run(graded).sam()
	g = graded.g
	byq = sam.by_qname()
	for q in ("q", "qr"):
		recs = byq[q]
		if outs is None:
			assert recs[0].line == base.by_qname()[q][0].line
		else:
			assert recs[0].line.split("\t")[:11] == base.by_qname()[q][0].line.split("\t")[:11]
		sec = recs[1:]
		assert all(r.is_secondary for r in sec)
		passing = g.secondaries_at(float(outs) if outs else 0.8)
		want = passing[:outn]
		assert len(sec) == len(want)
		for r, c in zip(sec, want):
			strand = c["strand"] if q == "q" else "+-"[c["strand"] == "+"]
			assert (r.rname, r.pos, r.strand, r.cigar, r.mapq, r.seq, r.qual) == (c["ctg"], c["pos"], strand, "150M", 0, "*", "*")
			assert r.flag == 0x100 | (0x10 if strand == "-" else 0)
			assert r.tags["NM"] == c["nm"] and r.tags["AS"] == r.tags["ms"] == c["score"]
	bsec = byq["b"][1:]
	assert len(bsec) == min(outn, 6)
	seen = set((r.rname, r.strand, r.pos) for r in byq["b"])
	assert seen <= set((c["ctg"], c["strand"], c["pos"]) for c in g.copies_b)
	assert len(seen) == len(byq["b"])


def test_outn_parses_suffixes(graded):
	"""--outn and --xa take numbers with k/m/g suffixes and decimals (kom_parse_num): --outn 0.001k is 1,
	--xa 0.006k is 6."""
	sam = _graded_run(graded, "--outn", "0.001k", "--xa", "0.006k").sam()
	byq = sam.by_qname()
	assert [len(byq[q]) for q in ("q", "qr", "b")] == [2, 2, 2]
	assert "XA" in byq["b"][0].tags


@pytest.mark.pinned
def test_pin_graded(graded):
	"""Pins the current output: SAM for the graded-repeat reads with --outn 10 --outs 0."""
	H.check_pinned("d_graded.sam", H.strip_pg(_graded_run(graded, "--outn", "10", "--outs", "0").stdout))
