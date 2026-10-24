"""Area A: reading queries and cutting them into batches.

bseq.c and the kseq parser are driven directly (tests/c/drv_a_input.c) and
through `minibwa map`: FASTA/FASTQ parsing, U-to-T conversion, comments and
`-y`, gzip and stdin, two-file and interleaved paired input, the grouping of
adjacent same-name reads into pairs, `-K` batching (checked against a model
of the documented rule using the per-batch "mapped N bp in M sequences"
lines on stderr), and the independence of single-end output from where
batches and sub-batches are cut.
"""

import gzip
import os

import pytest

import minibwa_helpers as H
import helpers_a as A


@pytest.fixture(scope="module")
def tiny(index_of, tmp_path_factory):
	return A.Tiny(index_of, str(tmp_path_factory.mktemp("tiny")))


def _records(text):
	"""Batches of (l_seq, name, seq, qual, comment) from drv_a_input's B/R lines."""
	batches = []
	for line in text.split("\n"):
		if line.startswith("B\t"):
			batches.append([])
		elif line.startswith("R\t"):
			f = line.split("\t", 5)
			batches[-1].append((int(f[1]), f[2], f[3], None if f[4] == "*" else f[4], None if f[5] == "*" else f[5]))
	return batches


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

# (raw text, name, sequence after U->T, quality, comment)
_PARSE_CASES = [
	(">a1 the first comment\nACGTU\nacgtu\nNNRY\n", "a1", "ACGTTacgttNNRY", None, "the first comment"),
	(">a2\tcomment after a tab\n\nGATTACA\n\n\nGG\n", "a2", "GATTACAGG", None, "comment after a tab"),
	(">empty_seq\n", "empty_seq", "", None, None),
	("@q1 x=1 y=2\nACGTN\n+\n!#%&I\n", "q1", "ACGTN", "!#%&I", "x=1 y=2"),
	("@q2\nUUUU\n+q2\nABCD\n", "q2", "TTTT", "ABCD", None),
	(">a3  two spaces\nTTTTGGGGCCCC\n", "a3", "TTTTGGGGCCCC", None, " two spaces"),
	("@q3 crlf\r\nACGT\r\n+\r\nIIII\r\n", "q3", "ACGT", "IIII", "crlf"),
	(">a4 crlf\r\nACGTAC\r\nGT\r\n", "a4", "ACGTACGT", None, "crlf"),
	("@q4\nAC\nGT\n+\nII\nII\n", "q4", "ACGT", "IIII", None),
	(">a5\n" + "ACGT" * 300 + "\n", "a5", "ACGT" * 300, None, None),
]


def _parse_input(tmp_path, gz=False):
	text = "".join(c[0] for c in _PARSE_CASES)
	path = str(tmp_path / ("in.fx.gz" if gz else "in.fx"))
	with (gzip.open(path, "wb") if gz else open(path, "wb")) as f:
		f.write(text.encode())
	return path


@pytest.mark.parametrize("with_qual,with_comment", [(1, 1), (0, 0), (1, 0)])
@pytest.mark.parametrize("gz", [False, True])
def test_bseq_parses_fasta_and_fastq(tmp_path, with_qual, with_comment, gz):
	"""mb_bseq_open/mb_bseq_read on mixed FASTA and FASTQ, plain or gzip: multi-line and
	empty-line FASTA, multi-line FASTQ, CRLF line ends, an empty sequence; U/u become
	T/t; quality and comment are kept only when asked for; names stop at the first
	whitespace and the comment is the rest of the line."""
	path = _parse_input(tmp_path, gz)
	r = H.assert_c_case("drv_a_input", "read", path, 100000000, with_qual, with_comment, 0, 40000, 1000000000)
	got = _records(r.stdout)
	assert len(got) == 1
	want = [(len(s), n, s, q if with_qual else None, c if with_comment else None) for _, n, s, q, c in _PARSE_CASES]
	assert got[0] == want


def test_bseq_reads_stdin_and_reports_open_failure(tmp_path):
	"""'-' reads standard input; a missing file cannot be opened."""
	path = _parse_input(tmp_path)
	with open(path, "rb") as f:
		data = f.read()
	r = H.run_c_case("drv_a_input", "read", ["-", 100000000, 1, 1, 0, 40000, 1000000000], stdin=data)
	assert r.passed, r.report()
	assert len(_records(r.stdout)[0]) == len(_PARSE_CASES)
	r = H.assert_c_case("drv_a_input", "open", path, str(tmp_path / "missing.fq"), "-")
	assert A.parse_rows(r.stdout, "O") == [[path, "1"], [str(tmp_path / "missing.fq"), "0"], ["-", "1"]]


def _qname_len(s):
	return len(s) - 2 if len(s) >= 3 and s[-1].isdigit() and s[-2] == "/" else len(s)


def test_qname_matching():
	"""mb_qname_len drops a trailing /<digit> from names of at least 3 characters;
	mb_qname_same compares names without it."""
	names = ["r1/1", "r1/2", "r1", "r1/12", "/1", "a/9", "a/x", "r1/", "ab", "ab/1/2", "a"]
	r = H.assert_c_case("drv_a_input", "qname", *names)
	for name, ln in A.parse_rows(r.stdout, "Q"):
		assert int(ln) == _qname_len(name), name
	n_same = 0
	for a, b, same in A.parse_rows(r.stdout, "S"):
		want = a[:_qname_len(a)] == b[:_qname_len(b)]
		assert int(same) == want, (a, b)
		n_same += want
	assert n_same > len(names)


# ---------------------------------------------------------------------------
# Batching
# ---------------------------------------------------------------------------


def bseq_batches(lens, names, chunk, min_cnt, max_chunk, frag):
	"""Record indices of the batches mb_bseq_read returns, per the rule in bseq.c: a batch
	ends once it holds at least max_chunk bases, or at least chunk bases and min_cnt
	records (max_chunk is raised to chunk; chunk <= 0 ends every batch after one record);
	in frag mode the following records with the same name (ignoring /1 /2) join the batch.
	A record read ahead to look for a mate starts the next batch."""
	max_chunk = max(max_chunk, chunk)
	out, i, n, carry = [], 0, len(lens), None
	while True:
		batch = [] if carry is None else [carry]
		carry = None
		while i < n:
			batch.append(i)
			i += 1
			size = sum(lens[j] for j in batch)
			if chunk <= 0 or max_chunk <= 0 or size >= max_chunk or (size >= chunk and len(batch) >= min_cnt):
				if frag and lens[batch[-1]] < 1000000:
					while i < n:
						j = i
						i += 1
						if names[j][:_qname_len(names[j])] == names[batch[-1]][:_qname_len(names[batch[-1]])]:
							batch.append(j)
						else:
							carry = j
							break
				break
		if not batch:
			return out
		out.append(batch)


def _batch_reads(seed, n, pairs_every=0, fixed=None):
	"""n records of 1-300 bp (or all `fixed` bp, so that batch sizes hit -K exactly);
	every pairs_every-th name is a pair with /1 and /2."""
	rng = H.Rng(seed)
	recs = []
	i = 0
	while len(recs) < n:
		ln = fixed or rng.randint(1, 300)
		if pairs_every and i % pairs_every == 0:
			recs.append(("p%d/1" % i, rng.dna(ln)))
			recs.append(("p%d/2" % i, rng.dna(fixed or rng.randint(1, 300))))
		else:
			recs.append(("s%d" % i, rng.dna(ln)))
		i += 1
	return recs


@pytest.mark.parametrize("chunk,min_cnt,max_chunk,frag,fixed", [
	(1000, 1, 1000, 0, None), (1000, 5, 3000, 0, None), (1000, 40000, 2500, 1, None), (2500, 3, 1000, 1, None),
	(0, 1, 0, 0, None), (0, 1, 0, 1, None), (1, 1, 1, 1, None), (100000, 1, 100000, 1, None),
	(300, 4, 1000, 0, 100), (300, 4, 1000, 1, 100), (400, 1, 400, 1, 100),
])
def test_bseq_batches(tmp_path, chunk, min_cnt, max_chunk, frag, fixed):
	"""mb_bseq_read cuts batches per the rule modelled in bseq_batches(), on single and
	interleaved paired records, also when a batch reaches -K exactly; no batch is empty
	and every record appears once, in order."""
	recs = _batch_reads(7, 120, pairs_every=3, fixed=fixed)
	path = A.write_reads_fasta(str(tmp_path / "b.fa"), recs)
	r = H.assert_c_case("drv_a_input", "read", path, chunk, 0, 0, frag, min_cnt, max_chunk)
	got = _records(r.stdout)
	want = bseq_batches([len(s) for _, s in recs], [n for n, _ in recs], chunk, min_cnt, max_chunk, frag)
	assert [[x[1] for x in b] for b in got] == [[recs[j][0] for j in b] for b in want]
	assert [x[2] for b in got for x in b] == [s for _, s in recs]


@pytest.mark.parametrize("chunk,fixed", [(0, None), (1, None), (700, None), (5000, None), (10 ** 9, None), (400, 100)])
def test_bseq_read_frag(tmp_path, chunk, fixed):
	"""mb_bseq_read_frag reads two files in lockstep, record i of each file adjacent, and ends
	a batch once it holds at least chunk bases (exactly chunk included); when one file runs
	out, the extra records of the other are dropped and reading stops."""
	r1 = _batch_reads(8, 40, fixed=fixed)
	r2 = _batch_reads(9, 46, fixed=fixed)
	f1 = A.write_reads_fasta(str(tmp_path / "1.fa"), r1)
	f2 = A.write_reads_fasta(str(tmp_path / "2.fa"), r2)
	r = H.assert_c_case("drv_a_input", "frag", chunk, 0, 0, f1, f2)
	got = _records(r.stdout)
	flat = [x[1] for b in got for x in b]
	assert flat == [n for p in zip(r1, r2) for n in (p[0][0], p[1][0])]
	want, size, cur = [], 0, []
	for a, b in zip(r1, r2):
		cur.append(a[0])
		cur.append(b[0])
		size += len(a[1]) + len(b[1])
		if size >= chunk:
			want.append(cur)
			cur, size = [], 0
	if cur:
		want.append(cur)
	assert [[x[1] for x in b] for b in got] == want


def _batch_stats(stderr):
	return A.batch_lines(stderr)


def _expect_stats(recs, batches):
	return [(sum(len(recs[j][1]) for j in b), len(b)) for b in batches]


def test_map_K_batches_single_file(tiny, tmp_path):
	"""minibwa map -K NUM1[,NUM2]: with fewer than 40,000 reads in a batch, a batch ends at
	NUM2 bases (NUM2 defaults to NUM1 and is raised to NUM1 when smaller); in the default
	paired-end mode the mate of a batch's last read joins it; with --pe=no it does not.
	NUM takes k/m/g suffixes. The batch sizes come from the per-batch stderr lines."""
	recs = _batch_reads(10, 90, pairs_every=4)
	fa = A.write_reads_fasta(str(tmp_path / "r.fa"), recs)
	lens, names = [len(s) for _, s in recs], [n for n, _ in recs]
	cases = [
		(["-K", "1000"], 1000, 1000, 1),
		(["-K", "1k"], 1000, 1000, 1),
		(["-K", "1000,3000"], 1000, 3000, 1),
		(["-K", "3000,1000"], 3000, 3000, 1),
		(["-K", "0.7k", "--pe=no"], 700, 700, 0),
		(["-K", "1"], 1, 1, 1),
		([], 100000000, 1000000000, 1),
	]
	for args, chunk, max_chunk, frag in cases:
		r = H.run_map(tiny.index, fa, args)
		want = bseq_batches(lens, names, chunk, 40000, max_chunk, frag)
		assert _batch_stats(r.stderr) == _expect_stats(recs, want), args


def test_map_K_batches_two_files(tiny, tmp_path):
	"""With two read files, a batch ends once its pairs hold at least NUM1 bases."""
	pairs = [("f" + n, x) for n, x in _batch_reads(12, 10, fixed=100)] + _batch_reads(11, 20)
	r1 = [(n, s) for n, s in pairs]
	r2 = [(n, H.revcomp(s)) for n, s in pairs]
	f1 = A.write_reads_fasta(str(tmp_path / "1.fa"), r1)
	f2 = A.write_reads_fasta(str(tmp_path / "2.fa"), r2)
	for k in (1, 800, 2500, 400):
		r = H.run_map(tiny.index, [f1, f2], ["-K", k])
		want, cur = [], [0, 0]
		for n, s in pairs:
			cur = [cur[0] + 2 * len(s), cur[1] + 2]
			if cur[0] >= k:
				want.append(tuple(cur))
				cur = [0, 0]
		if cur[1]:
			want.append(tuple(cur))
		assert _batch_stats(r.stderr) == want, k


def test_map_batch_needs_40000_reads(tiny, tmp_path):
	"""A batch that reaches NUM1 bases still grows to 40,000 reads unless it reaches NUM2."""
	path = str(tmp_path / "many.fa")
	with open(path, "w") as f:
		for i in range(40007):
			f.write(">r%d\nACG\n" % i)
	r = H.run_map(tiny.index, path, ["-K", "1000,1g", "-f"])
	assert _batch_stats(r.stderr) == [(120000, 40000), (21, 7)]
	r = H.run_map(tiny.index, path, ["-K", "1000,30000", "-f"])
	assert _batch_stats(r.stderr)[:2] == [(30000, 10000), (30000, 10000)]


def test_batching_does_not_change_single_end_output(tiny, tmp_path):
	"""Single-end output is the same whether every read is its own batch (-K 1) or all
	reads share one batch cut into sub-batches of 24 reads, at any -t."""
	recs = [(n, s) for n, s in tiny.reads]
	fq = H.write_fastq(str(tmp_path / "r.fq"), [(n, s, "I" * len(s)) for n, s in recs])
	whole = H.run_map(tiny.index, fq, ["-k", "15"])
	assert len(_batch_stats(whole.stderr)) == 1
	assert len(H.parse_sam(whole.stdout).primaries()) == len(recs)
	one = H.run_map(tiny.index, fq, ["-k", "15", "-K", "1", "--pe=no"], threads=1)
	assert len(_batch_stats(one.stderr)) == len(recs)
	assert H.sam_body(one.stdout) == H.sam_body(whole.stdout)
	two = H.run_map(tiny.index, fq, ["-k", "15", "-K", "1"])
	assert len(_batch_stats(two.stderr)) == (len(recs) + 2) // 2  # a read is read ahead to look for its mate
	assert H.sam_body(two.stdout) == H.sam_body(whole.stdout)


def test_subbatch_limits_do_not_change_hits(tiny):
	"""mb_map_batch gives every read the same hits for any sub-batch limits (sb_seq 1-1000
	reads, sb_len 1 bp-1 Mbp) as with the defaults and as mb_map one read at a time."""
	H.assert_c_case("drv_a_input", "subbatch", tiny.index, tiny.reads_fa)


def test_option_defaults():
	"""mb_opt_init: -k 19, -c 250 and -K 100m,1g as in the usage text, paired-end mode on;
	the lr preset batches 500m bases."""
	H.assert_c_case("drv_a_input", "opt")


# ---------------------------------------------------------------------------
# minibwa map input forms
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def pe(tiny, tmp_path_factory):
	"""60 FR pairs from the tiny genome's long contigs and the SAM of their two-file run."""
	g = H.Genome(tiny.contigs[:2])
	pairs = H.simulate_pairs(g, 60, read_length=100, insert_mean=250, insert_sd=20, seed=31, sub_rate=0.005,
		indel_rate=0, contigs=["t1", "t2"])
	d = str(tmp_path_factory.mktemp("pe"))
	f1, f2 = H.write_pairs(os.path.join(d, "r1.fq"), os.path.join(d, "r2.fq"), pairs)
	sam = H.run_map(tiny.index, [f1, f2], ["-I", "250,20"])
	return pairs, d, f1, f2, sam


def test_input_forms_give_the_same_output(tiny, pe):
	"""The same pairs as two FASTQ files, two gzip files, one interleaved file with or without
	/1 /2 suffixes, and interleaved on stdin give the same SAM records; two FASTA files
	give the same records with QUAL '*'."""
	pairs, d, f1, f2, sam = pe
	body = H.sam_body(sam.stdout)
	assert len(H.parse_sam(sam.stdout).records) >= 120
	g1, g2 = H.write_pairs(os.path.join(d, "r1.fq.gz"), os.path.join(d, "r2.fq.gz"), pairs)
	assert H.sam_body(H.run_map(tiny.index, [g1, g2], ["-I", "250,20"]).stdout) == body
	il = H.write_interleaved(os.path.join(d, "il.fq"), pairs, suffix=True)
	assert H.sam_body(H.run_map(tiny.index, il, ["-I", "250,20"]).stdout) == body
	il2 = H.write_interleaved(os.path.join(d, "il2.fq.gz"), pairs)
	assert H.sam_body(H.run_map(tiny.index, il2, ["-I", "250,20"]).stdout) == body
	with open(il, "rb") as f:
		data = f.read()
	assert H.sam_body(H.run_map(tiny.index, "-", ["-I", "250,20"], stdin=data).stdout) == body
	a1, a2 = H.write_pairs(os.path.join(d, "r1.fa"), os.path.join(d, "r2.fa"), pairs)
	fa = H.parse_sam(H.run_map(tiny.index, [a1, a2], ["-I", "250,20"]).stdout)
	fq = H.parse_sam(sam.stdout)
	assert len(fa.records) == len(fq.records)
	for x, y in zip(fa.records, fq.records):
		assert x.qual == "*" and x.line.split("\t")[:10] == y.line.split("\t")[:10]


def test_read_grouping_and_name_trimming(tiny, tmp_path):
	"""In one input file, adjacent reads with the same name (ignoring a /<digit> suffix) form
	a pair (flags 0x1, 0x40/0x80) and lose differing /1 /2 suffixes in QNAME; other reads
	stay single with their names as given; mixed single and paired reads work; with
	--pe=no nothing is paired."""
	g = dict(tiny.contigs)
	t1 = g["t1"]
	recs = [
		("pa/1", t1[100:200]), ("pa/2", H.revcomp(t1[300:400])),
		("single/1", t1[500:600]),
		("pb", t1[700:800]), ("pb", H.revcomp(t1[900:1000])),
		("pc/1", t1[1000:1100]), ("pc/1", H.revcomp(t1[1150:1250])),
		("lone", t1[1300:1400]),
		("pd/2", t1[1600:1700]), ("pd/1", H.revcomp(t1[1750:1850])),
	]
	fa = A.write_reads_fasta(str(tmp_path / "g.fa"), recs)
	sam = H.run_map(tiny.index, fa).sam()
	got = [(x.qname, x.flag & 0xC1) for x in sam.records if x.is_primary]
	assert got == [("pa", 0x41), ("pa", 0x81), ("single/1", 0), ("pb", 0x41), ("pb", 0x81),
		("pc/1", 0x41), ("pc/1", 0x81), ("lone", 0), ("pd", 0x41), ("pd", 0x81)]
	assert H.pair_invariant_errors(sam) == []
	sam = H.run_map(tiny.index, fa, ["--pe=no"]).sam()
	assert [(x.qname, x.flag & 0xC1) for x in sam.records if x.is_primary] == [(n, 0) for n, _ in recs]


def test_two_files_pair_by_position_and_name(tiny, tmp_path):
	"""Two files are read in lockstep; records at the same index pair up when their names
	match (ignoring /1 /2) and stay single otherwise; extra records in the longer file are
	dropped with a warning."""
	t1 = dict(tiny.contigs)["t1"]
	r1 = [("x/1", t1[100:200]), ("y", t1[400:500]), ("z", t1[700:800]), ("extra1", t1[1000:1100])]
	r2 = [("x/2", H.revcomp(t1[250:350])), ("w", H.revcomp(t1[550:650])), ("z", H.revcomp(t1[850:950]))]
	f1 = A.write_reads_fasta(str(tmp_path / "1.fa"), r1)
	f2 = A.write_reads_fasta(str(tmp_path / "2.fa"), r2)
	res = H.run_map(tiny.index, [f1, f2])
	got = [(x.qname, x.flag & 0xC1) for x in res.sam().records if x.is_primary]
	assert got == [("x", 0x41), ("x", 0x81), ("y", 0), ("w", 0), ("z", 0x41), ("z", 0x81)]
	assert "different number of records" in res.stderr


def test_copy_comment(tiny, tmp_path):
	"""-y appends each read's FASTA/FASTQ comment verbatim as the last field of its SAM
	records (mapped and unmapped) and of the PAF lines of mapped reads; without -y no
	comment is written. (Base minibwa writes no comment on the PAF line of an unmapped
	read, so that line is not checked.)"""
	t1 = dict(tiny.contigs)["t1"]
	recs = [("c1", t1[100:220], "I" * 120, "BC:Z:ACGT\tXX:i:7"), ("c2", H.revcomp(t1[600:700]), "I" * 100, "plain words"),
		("c3", "ACGT" * 20, "I" * 80, "unmapped one"), ("c4", t1[1000:1100], "I" * 100, None)]
	fq = H.write_fastq(str(tmp_path / "c.fq"), recs)
	want = dict((r[0], r[3]) for r in recs)
	for fmt in ([], ["-f"]):
		with_c = H.run_map(tiny.index, fq, fmt + ["-y", "--pe=no"]).stdout
		without = H.run_map(tiny.index, fq, fmt + ["--pe=no"]).stdout
		lines = [l for l in with_c.split("\n") if l and not l.startswith("@")]
		plain = [l for l in without.split("\n") if l and not l.startswith("@")]
		assert len(lines) == len(plain) >= 4
		n_unmapped = 0
		for a, b in zip(lines, plain):
			name = a.split("\t")[0]
			c = want[name]
			if fmt and b.split("\t")[5] == "*":
				n_unmapped += 1
				assert a.startswith(b), name
				continue
			assert a == (b + "\t" + c if c else b), name
		assert n_unmapped == (1 if fmt else 0)
	# a comment holding SAM tags becomes tags of the record
	out = H.run_map(tiny.index, fq, ["-y", "--pe=no", "-x", "sr"]).stdout
	rec = H.SamRecord([l for l in out.split("\n") if l.startswith("c1\t")][0])
	assert rec.tags["BC"] == "ACGT" and rec.tags["XX"] == 7


def test_parse_warnings(tiny, tmp_path):
	"""A FASTQ record whose quality is shorter than its sequence stops reading with a
	warning naming the record before it (or the first record); reads before it are mapped.
	An empty read name draws a warning."""
	t1 = dict(tiny.contigs)["t1"]
	good = "".join("@g%d\n%s\n+\n%s\n" % (i, t1[100 * i:100 * i + 80], "I" * 80) for i in range(1, 4))
	path = str(tmp_path / "bad.fq")
	with open(path, "w") as f:
		f.write(good + "@bad\nACGTACGT\n+\nIII\n")
	r = H.run_map(tiny.index, path, ["--pe=no"])
	assert "failed to parse the FASTA/FASTQ record next to 'g3'" in r.stderr
	assert [x.qname for x in r.sam().records if x.is_primary] == ["g1", "g2", "g3"]
	with open(path, "w") as f:
		f.write("@bad\nACGTACGT\n+\nIII\n")
	r = H.run_map(tiny.index, path)
	assert "failed to parse the first FASTA/FASTQ record" in r.stderr
	with open(path, "w") as f:
		f.write(">\n%s\n>named\n%s\n" % (t1[100:200], t1[300:400]))
	r = H.run_map(tiny.index, path, ["-f"])
	assert "empty sequence name" in r.stderr
	assert [l.split("\t")[0] for l in r.stdout.split("\n") if l][-1] == "named"
