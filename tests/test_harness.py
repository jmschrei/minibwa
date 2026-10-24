"""Tests of the test harness itself: random numbers, genomes, the read
simulator's truth, tag recomputation, the parsers and the out-of-tree build."""

import gzip
import os
import shutil

import pytest

import conftest
import minibwa_build
import minibwa_helpers as H


# --- random numbers and genomes -------------------------------------------

def test_rng_matches_reference_splitmix64():
	# Reference outputs of splitmix64 seeded with 0 (Vigna's splitmix64.c).
	r = H.Rng(0)
	assert [r.u64() for _ in range(3)] == [0xE220A8397B1DCDAF, 0x6E789E6AA1B965F4, 0x06C45D188009454F]


def test_rng_helpers_stay_in_range():
	r = H.Rng(5)
	assert all(0 <= r.below(7) < 7 for _ in range(1000))
	assert all(3 <= r.randint(3, 5) <= 5 for _ in range(1000))
	assert all(0.0 <= r.random() < 1.0 for _ in range(1000))
	x = [r.normal_int(300, 30) for _ in range(4000)]
	mean = sum(x) / float(len(x))
	sd = (sum((v - mean) ** 2 for v in x) / len(x)) ** 0.5
	assert abs(mean - 300) < 3 and abs(sd - 30) < 3


def test_synth_genome_is_stable():
	"""The standard genome must not change across Python versions or platforms: pins depend on it."""
	g = conftest.make_synth()
	assert g.names == ["chr1", "chr2", "short1", "short2"]
	assert [g.length(n) for n in g.names] == [60000, 40000, 72, 82]
	assert g.digest() == "f74c1a0bc11acba678bd2e2af932d72a73b911b6b7f19c6810803f1d91628364"


def test_genome_features_are_what_they_claim():
	g = conftest.make_synth()
	kinds = {}
	for f in g.features:
		kinds[f.kind] = kinds.get(f.kind, 0) + 1
		s = g[f.contig][f.start:f.end]
		if f.kind == "N":
			assert set(s) == {"N"}
		elif f.kind == "lowcomp":
			assert s == (f.info * len(s))[:len(s)]
	assert kinds == {"repeat": 6, "lowcomp": 3, "N": 2}
	fams = {}
	for f in g.features:
		if f.kind == "repeat":
			s = g[f.contig][f.start:f.end]
			fams.setdefault(f.info[0], []).append(s if f.info[2] == "+" else H.revcomp(s))
	for copies in fams.values():  # divergence 0: identical copies
		assert len(set(copies)) == 1
	assert not g.is_unique(g.features[0].contig, g.features[0].start, g.features[0].start + 1)


# --- read simulator ------------------------------------------------------

def _apply_truth(ref, seg, fwd):
	"""Rebuild the read (reference orientation) from the reference, CIGAR and edits."""
	subs = dict((e[0], e[3]) for e in seg.edits if e[1] == "X")
	ins = dict((e[0], e[3]) for e in seg.edits if e[1] == "I")
	out, r = [], seg.pos
	for n, op in H.parse_cigar(seg.cigar):
		if op == "M":
			for _ in range(n):
				out.append(subs.get(r, ref[r].upper() if ref[r].upper() in "ACGT" else "N"))
				r += 1
		elif op == "I":
			out.append(ins[r])
		elif op == "D":
			r += n
	assert r == seg.end
	return "".join(out)


def test_simulated_reads_match_their_truth():
	g = conftest.make_synth()
	reads = H.simulate_reads(g, 300, seed=2, sub_rate=0.02, indel_rate=0.002)
	assert len(set(r.name for r in reads)) == 300
	n_ind = 0
	for r in reads:
		assert len(r.seq) == len(r.qual) == 150
		assert H.cigar_query_span(r.cigar) == 150
		assert H.cigar_ref_span(r.cigar) == r.end - r.pos
		assert _apply_truth(g[r.contig], r.segments[0], r.fwd_seq) == r.fwd_seq
		assert r.seq == (r.fwd_seq if r.strand == "+" else H.revcomp(r.fwd_seq))
		# bases copied from an N run are not edits, but the SAM spec counts them in NM
		assert H.compute_nm(g[r.contig], r.pos, r.fwd_seq, r.cigar) <= r.nm + r.fwd_seq.count("N")
		n_ind += "I" in r.cigar or "D" in r.cigar
	assert n_ind > 10


def test_unique_reads_avoid_features():
	g = conftest.make_synth()
	for r in H.simulate_reads(g, 200, seed=3, unique=True):
		assert g.is_unique(r.contig, r.pos, r.end)


def test_simulated_pairs_are_fr_with_the_drawn_insert():
	g = conftest.make_synth()
	pairs = H.simulate_pairs(g, 200, seed=4, insert_mean=350, insert_sd=40)
	for p in pairs:
		assert p.r1.name == p.r2.name == p.name
		assert p.r1.contig == p.r2.contig == p.contig
		left, right = (p.r1, p.r2) if p.strand == "+" else (p.r2, p.r1)
		assert left.strand == "+" and right.strand == "-"
		assert left.pos == p.frag_start and right.end == p.frag_end
		for r in p:
			assert _apply_truth(g[r.contig], r.segments[0], r.fwd_seq) == r.fwd_seq
	ins = [p.insert for p in pairs]
	mean = sum(ins) / float(len(ins))
	assert 330 < mean < 370


def test_hic_pairs_and_chimeras():
	g = conftest.make_synth()
	pairs = H.simulate_hic_pairs(g, 100, seed=5, chimeric_frac=0.5, min_distance=20000)
	n_chim = 0
	for p in pairs:
		a, b = p.r1.segments[0], p.r2.segments[0]
		assert a.contig != b.contig or abs(a.pos - b.pos) >= 20000
		if p.r1.is_chimeric:
			n_chim += 1
			s1, s2 = p.r1.segments
			assert s1.qstart == 0 and s1.qend == s2.qstart and s2.qend == len(p.r1.seq)
			for s in (s1, s2):
				part = p.r1.seq[s.qstart:s.qend]
				fwd = part if s.strand == "+" else H.revcomp(part)
				assert _apply_truth(g[s.contig], s, fwd) == fwd
	assert 30 < n_chim < 70


def test_long_reads_have_variable_length():
	g = conftest.make_synth()
	reads = H.simulate_long_reads(g, 20, seed=6, length=(1000, 3000))
	assert len(set(len(r) for r in reads)) > 10
	for r in reads:
		assert 1000 <= len(r) <= 3000
		assert _apply_truth(g[r.contig], r.segments[0], r.fwd_seq) == r.fwd_seq


# --- FASTA/FASTQ writers and readers --------------------------------------

def test_fastx_round_trip_plain_and_gzip(tmp_path):
	g = conftest.make_synth()
	reads = H.simulate_reads(g, 10, seed=7)
	for name in ("r.fq", "r.fq.gz", "r.fa", "r.fa.gz"):
		p = H.write_fastx(str(tmp_path / name), reads)
		back = H.read_fastx(p)
		assert [x[0] for x in back] == [r.name for r in reads]
		assert [x[1] for x in back] == [r.seq for r in reads]
		if ".fq" in name:
			assert [x[2] for x in back] == [r.qual for r in reads]
	with open(str(tmp_path / "r.fq.gz"), "rb") as f:
		assert f.read(2) == b"\x1f\x8b"
	pairs = H.simulate_pairs(g, 5, seed=8)
	H.write_interleaved(str(tmp_path / "i.fq"), pairs, suffix=True)
	names = [x[0] for x in H.read_fastx(str(tmp_path / "i.fq"))]
	assert names == [n for p in pairs for n in (p.name + "/1", p.name + "/2")]


def test_gzip_output_is_byte_stable(tmp_path):
	recs = [("a", "ACGT", "IIII")]
	H.write_fastq(str(tmp_path / "x.fq.gz"), recs)
	H.write_fastq(str(tmp_path / "y.fq.gz"), recs)
	with open(str(tmp_path / "x.fq.gz"), "rb") as f, open(str(tmp_path / "y.fq.gz"), "rb") as g:
		assert f.read() == g.read()


# --- tag recomputation, worked by hand ------------------------------------

REF = "ACGTACGTAC"


@pytest.mark.parametrize("seq,cigar,nm,md,cs", [
	("ACGTACGTAC", "10M", 0, "10", ":10"),
	("ACGAACGTAC", "10M", 1, "3T6", ":3*ta:6"),
	("TCGTACGTAC", "10M", 1, "0A9", "*at:9"),
	("ACGTACGTAG", "10M", 1, "9C0", ":9*cg"),
	("ACGCGTAC", "3M2D5M", 2, "3^TA5", ":3-ta:5"),
	("ACGCGTAA", "3M2D5M", 3, "3^TA4C0", ":3-ta:4*ca"),
	("ACGTTTACGTAC", "4M2I6M", 2, "10", ":4+tt:6"),
	("GGACGTACGTAC", "2S10M", 0, "10", ":10"),
	("ACGTACGTNC", "10M", 1, "8A1", ":8*an:1"),
	("ACGAGTAC", "3M2D5M", 3, "3^TA0C4", ":3-ta*ca:4"),
	("ACGTTCGTAC", "4=1X5=", 1, "4A5", ":4*at:5"),
])
def test_tag_recomputation_hand_worked(seq, cigar, nm, md, cs):
	assert H.compute_nm(REF, 0, seq, cigar) == nm
	assert H.compute_md(REF, 0, seq, cigar) == md
	assert H.compute_cs(REF, 0, seq, cigar) == cs


def test_tag_recomputation_offsets_and_meth():
	ref = "TTTTACGTACGTAC"
	assert H.compute_nm(ref, 4, "ACGTACGTAC", "10M") == 0
	assert H.compute_md(ref, 4, "ACGTATGTAC", "10M") == "5C4"
	assert H.compute_nm(ref, 4, "ACGTATGTAC", "10M", meth="C2T") == 0
	assert H.compute_nm(ref, 4, "ACGTATGTAC", "10M", meth="G2A") == 1
	assert H.compute_cs(ref, 4, "ACGTACGTAC", "10M", long_form=True) == "=ACGTACGTAC"
	assert H.cigar_ref_span("5S10M2D3M1I4M3H") == 19
	assert H.cigar_query_span("5S10M2D3M1I4M3H") == 23
	assert H.cigar_query_span("5S10M2D3M1I4M3H", include_hard=True) == 26
	assert H.md_from_minibwa("149C0") == "149C" and H.md_from_minibwa("3^TA0") == "3^TA"
	assert H.md_from_minibwa("10") == "10" and H.md_from_minibwa("5C0A0") == "5C0A"


def test_n_never_matches_in_nm_and_md():
	assert H.compute_nm("ACNT", 0, "ACNT", "4M") == 1
	assert H.compute_md("ACNT", 0, "ACNT", "4M") == "2N1"
	assert H.compute_md("ACNT", 0, "ACNT", "4M", ambiguous_equal=True) == "4"
	assert H.compute_cs("ACNT", 0, "ACNT", "4M") == ":4"


# --- parsers -------------------------------------------------------------

SAM = (
	"@HD\tVN:1.6\tSO:unsorted\n"
	"@SQ\tSN:c1\tLN:100\n"
	"@RG\tID:g1\tSM:s\n"
	"@PG\tID:minibwa\tPN:minibwa\tVN:x\tCL:minibwa map a b\n"
	"q1\t99\tc1\t5\t60\t2S8M\t=\t40\t50\tAACGTACGTA\tIIIIIIIIII\tNM:i:0\tAS:i:-3\tXf:f:1.5e2\tXa:A:P\tXh:H:1AFF\tXb:B:i,1,-2,3\tMC:Z:10M\n"
	"q1\t147\tc1\t40\t60\t10M\t=\t5\t-50\tACGTACGTAC\t*\n"
	"q2\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t*\n"
)


def test_sam_parser_types_and_flags():
	s = H.parse_sam(SAM)
	assert s.hd == {"VN": "1.6", "SO": "unsorted"}
	assert s.sq == [("c1", 100)]
	assert s.rg == [{"ID": "g1", "SM": "s"}]
	assert s.pg[0]["CL"] == "minibwa map a b"
	r = s.records[0]
	assert (r.qname, r.flag, r.rname, r.pos, r.mapq, r.cigar, r.rnext, r.pnext, r.tlen) == ("q1", 99, "c1", 5, 60, "2S8M", "=", 40, 50)
	assert r.tags == {"NM": 0, "AS": -3, "Xf": 150.0, "Xa": "P", "Xh": "1AFF", "Xb": ("i", [1, -2, 3]), "MC": "10M"}
	assert r.tag_types["Xf"] == "f" and r.tag_order[0] == "NM"
	assert r.is_paired and r.is_proper_pair and r.is_read1 and r.mate_reverse and not r.is_reverse and r.is_primary
	assert r.pos0 == 4 and r.ref_end0 == 12 and r.left_clip == 2 and r.unclipped_start0 == 2 and r.query_length == 10
	r2 = s.records[1]
	assert r2.is_reverse and r2.is_read2 and r2.segment == 2 and r2.strand == "-"
	assert s.records[2].is_unmapped and s.records[2].segment == 0
	assert set(s.primaries()) == {("q1", 1), ("q1", 2), ("q2", 0)}
	assert H.strip_pg(SAM).count("@PG") == 0 and H.sam_body(SAM).count("\n") == 3


@pytest.mark.parametrize("bad", [
	"q\t0\tc1\t1\t60\t4M\t*\t0\t0\tACGT\t*\tNM:i:0\tNM:i:1\n",  # duplicate tag
	"q\t0\tc1\t1\t60\t5M\t*\t0\t0\tACGT\t*\n",  # CIGAR and SEQ lengths differ
	"q\t0\tc1\t1\t60\t4M\t*\t0\t0\tACGT\tIII\n",  # QUAL length
	"q\t0\tc1\t1\t300\t4M\t*\t0\t0\tACGT\t*\n",  # MAPQ range
	"q\t0\tc1\t1\t60\t4M\t*\t0\t0\tACGT\t*\tNM:i:x\n",  # bad integer
	"q\t0\tc1\t1\t60\t4M\t*\t0\t0\tACGT\n",  # too few fields
	"q\t0\tc1\t1\t60\t4M\t*\t0\t0\tACGT\t*",  # no final newline
	"q\t0\tc1\t1\t60\t4Q\t*\t0\t0\tACGT\t*\n",  # bad CIGAR op
])
def test_sam_parser_rejects_malformed(bad):
	with pytest.raises(H.SamFormatError):
		H.parse_sam(bad)


def test_paf_parser():
	text = ("r1\t150\t0\t150\t+\tc1\t1000\t10\t160\t148\t150\t60\ttp:A:P\tcm:i:5\tcg:Z:150M\n"
		"r2\t100\t*\t*\t*\t*\t*\t*\t*\t0\t0\t0\n")
	p = H.parse_paf(text)
	assert (p[0].qname, p[0].qlen, p[0].qstart, p[0].qend, p[0].strand, p[0].tname, p[0].tlen, p[0].tstart, p[0].tend,
		p[0].nmatch, p[0].alen, p[0].mapq) == ("r1", 150, 0, 150, "+", "c1", 1000, 10, 160, 148, 150, 60)
	assert p[0].is_primary and p[0].tags["cg"] == "150M"
	assert p[1].is_unmapped and p[1].qlen == 100 and p[1].tstart is None


def test_pin_helper(tmp_path, monkeypatch):
	monkeypatch.setattr(H, "PINS_DIR", str(tmp_path))
	monkeypatch.setenv("MINIBWA_TEST_UPDATE_PINS", "1")
	H.check_pinned("a.txt", "x\ny\n")
	H.check_pinned("b.txt.gz", "x\ny\n")
	monkeypatch.delenv("MINIBWA_TEST_UPDATE_PINS")
	H.check_pinned("a.txt", "x\ny\n")
	H.check_pinned("b.txt.gz", "x\ny\n")
	with gzip.open(str(tmp_path / "b.txt.gz"), "rt") as f:
		assert f.read() == "x\ny\n"
	with pytest.raises(AssertionError, match="differs from pin"):
		H.check_pinned("a.txt", "x\nz\n")
	with pytest.raises(AssertionError, match="missing pin"):
		H.check_pinned("c.txt", "x\n")


# --- runner and build ----------------------------------------------------

def test_run_minibwa_reports_exit_status():
	r = H.run_minibwa(["no-such-command"])
	assert r.returncode != 0 and "unknown command" in r.stderr
	with pytest.raises(AssertionError, match="exited with"):
		r.check()
	assert H._has_thread_opt(["-t", "4"]) and H._has_thread_opt(["-t4"]) and not H._has_thread_opt(["-x", "sr"])


def test_source_hash_tracks_sources_only(tmp_path):
	src = minibwa_build.src_dir()
	dst = str(tmp_path / "src")
	os.makedirs(dst)
	for rel in minibwa_build.source_files(src):
		os.makedirs(os.path.dirname(os.path.join(dst, rel)), exist_ok=True)
		shutil.copy2(os.path.join(src, rel), os.path.join(dst, rel))
	h0 = minibwa_build.source_hash(dst)
	assert h0 == minibwa_build.source_hash(src)
	os.makedirs(os.path.join(dst, "tests"))
	with open(os.path.join(dst, "tests", "x.py"), "w") as f:
		f.write("# not a source file\n")
	assert minibwa_build.source_hash(dst) == h0
	with open(os.path.join(dst, "seed.c"), "a") as f:
		f.write("\n")
	assert minibwa_build.source_hash(dst) != h0


def test_make_args_change_the_build_key(monkeypatch):
	"""MINIBWA_TEST_MAKE_ARGS reaches the build key (so an ASan or coverage build gets its own cache entry);
	without it the key is unchanged."""
	monkeypatch.delenv("MINIBWA_TEST_MAKE_ARGS", raising=False)
	k0 = minibwa_build.build_key()
	monkeypatch.setenv("MINIBWA_TEST_MAKE_ARGS", "asan=1 'MALLOC_O='")
	assert minibwa_build.make_args() == ["asan=1", "MALLOC_O="]
	assert minibwa_build.build_key() != k0
	monkeypatch.setenv("MINIBWA_TEST_MAKE_ARGS", "")
	assert minibwa_build.build_key() == k0


def test_out_of_tree_build_leaves_tree_clean(minibwa_bin, tree_before):
	"""Building the binary, libminibwa.a and a C driver adds no build product to the source tree."""
	minibwa_build.c_driver("drv_pilot")
	src = minibwa_build.src_dir()
	if os.environ.get("MINIBWA_BIN") is None:
		assert not os.path.realpath(minibwa_bin).startswith(os.path.realpath(src) + os.sep)
	assert minibwa_build.tree_products(src) == tree_before
