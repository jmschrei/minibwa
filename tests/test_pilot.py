"""Pilot tests: whole-program sanity on the shipped chrM data, simulated
single-end and paired-end reads against their truth, and one libminibwa API
case through a C driver."""

import shutil
import subprocess

import pytest

import minibwa_helpers as H


def _chrm_reads(chrm):
	reads = {}
	for fn in (chrm.read1, chrm.read2):
		for name, seq, _, _ in H.read_fastx(fn):
			reads.setdefault(name, []).append(seq)
	return reads


@pytest.fixture(scope="module")
def chrm_pe(chrm):
	return H.run_map(chrm.index, [chrm.read1, chrm.read2])


def test_chrm_pe_sam_sanity(chrm, chrm_pe):
	"""Header from the genome; every pair present with one primary per mate; SEQ is the input read;
	CIGAR spans the read and stays on the contig; NM matches the reference; mate fields agree."""
	sam = chrm_pe.sam()
	assert sam.hd.get("VN") == "1.6"
	assert sam.sq == [(n, len(s)) for n, s in chrm.genome.items()] == [("chrM", 16569)]
	assert len(sam.pg) == 1 and sam.pg[0]["ID"] == "minibwa" and " map " in sam.pg[0]["CL"]
	reads = _chrm_reads(chrm)
	assert len(reads) == 1000 and all(len(v) == 2 for v in reads.values())
	assert H.sam_invariant_errors(sam, chrm.genome, reads) == []
	assert H.pair_invariant_errors(sam) == []
	prim = sam.primaries()
	assert set(prim) == set((n, s) for n in reads for s in (1, 2))
	# these are real human chrM reads: nearly all should pair properly
	n_proper = sum(r.is_proper_pair for r in prim.values())
	assert n_proper >= 0.95 * len(prim)


def test_chrm_tlen(chrm_pe):
	"""TLEN of every pair on one contig is the 5'-to-5' template length (the SAM spec leaves
	the definition to the implementation; minibwa uses its TLEN#1). For an FR pair whose
	forward mate starts and ends no later than the reverse mate, it also equals the
	leftmost-to-rightmost length (TLEN#2)."""
	prim = chrm_pe.sam().primaries()
	n_fr = 0
	for (q, s), r in prim.items():
		m = prim[(q, 3 - s)]
		if r.is_unmapped or m.is_unmapped or r.rname != m.rname:
			assert r.tlen == 0
			continue
		assert r.tlen == H.tlen_5prime(r, m), (r.line[:150], m.line[:150])
		f, v = (r, m) if not r.is_reverse else (m, r)
		if r.is_reverse != m.is_reverse and f.pos0 <= v.pos0 and f.ref_end0 <= v.ref_end0:
			assert r.tlen == H.tlen_outer(r, m), (r.line[:150], m.line[:150])
			n_fr += 1
	assert n_fr > 1800


def test_chrm_paf_sanity(chrm):
	"""PAF columns are in range, cg:Z spans the stated query and target intervals, NM matches."""
	res = H.run_map(chrm.index, [chrm.read1, chrm.read2], ["-f"])
	paf = res.paf()
	reads = _chrm_reads(chrm)
	seen = set()
	for p in paf:
		name, seg = p.qname.rsplit("/", 1)
		read = reads[name][int(seg) - 1]
		seen.add(p.qname)
		assert p.qlen == len(read)
		if p.is_unmapped:
			continue
		assert (p.tname, p.tlen) == ("chrM", 16569)
		assert 0 <= p.qstart < p.qend <= p.qlen and 0 <= p.tstart < p.tend <= p.tlen
		assert 0 <= p.nmatch <= p.alen and 0 <= p.mapq <= 60
		cg = p.tags["cg"]
		assert H.cigar_query_span(cg) == p.qend - p.qstart and H.cigar_ref_span(cg) == p.tend - p.tstart
		s = read if p.strand == "+" else H.revcomp(read)
		qs = p.qstart if p.strand == "+" else p.qlen - p.qend
		assert p.tags["NM"] == H.compute_nm(chrm.genome["chrM"], p.tstart, s[qs:qs + p.qend - p.qstart], cg)
	assert seen == set("%s/%d" % (n, i) for n in reads for i in (1, 2))


@pytest.fixture(scope="module")
def se_run(synth, tmp_path_factory):
	reads = H.simulate_reads(synth.genome, 300, seed=21, unique=True)
	fq = H.write_fastq(str(tmp_path_factory.mktemp("se") / "se.fq"), reads)
	return reads, H.run_map(synth.index, fq, ["-b", "MD"])


def test_synth_se_unique_reads_land_at_truth(synth, se_run):
	"""Reads from unique sequence: true contig and strand, unclipped start within 3 bp of the truth,
	MAPQ >= 20; CIGAR, SEQ and NM consistent with the genome; MD equals its recomputation."""
	reads, res = se_run
	sam = res.sam()
	truth = dict((r.name, r) for r in reads)
	assert H.sam_invariant_errors(sam, synth.genome, dict((r.name, [r.seq]) for r in reads)) == []
	prim = sam.primaries()
	assert len(prim) == len(reads)
	for (q, _), rec in prim.items():
		err = H.placement_error(rec, truth[q])
		assert err is not None and err <= 3, (rec.line[:150], truth[q])
		assert rec.mapq >= 20
		assert H.cigar_query_span(rec.cigar) == len(truth[q].seq)
		md = H.compute_md(synth.genome[rec.rname], rec.pos0, rec.seq, rec.cigar)
		# base minibwa omits the final 0 after a terminal mismatch, which the SAMtags MD regex requires
		assert rec.tags["MD"] in (md, H.md_from_minibwa(md))


def test_synth_pe_unique_pairs_land_at_truth(synth, tmp_path):
	"""FR pairs from unique sequence: both mates at the truth, flag 0x2 on both, TLEN equal to
	the 5'-to-5' length and within 10 bp of the simulated insert."""
	pairs = H.simulate_pairs(synth.genome, 300, seed=22, insert_mean=300, insert_sd=30, unique=True)
	f1, f2 = H.write_pairs(str(tmp_path / "r1.fq"), str(tmp_path / "r2.fq"), pairs)
	sam = H.run_map(synth.index, [f1, f2]).sam()
	assert H.sam_invariant_errors(sam, synth.genome, dict((p.name, [p.r1.seq, p.r2.seq]) for p in pairs)) == []
	assert H.pair_invariant_errors(sam) == []
	prim = sam.primaries()
	for p in pairs:
		a, b = prim[(p.name, 1)], prim[(p.name, 2)]
		for rec, read in ((a, p.r1), (b, p.r2)):
			err = H.placement_error(rec, read)
			assert err is not None and err <= 3, (rec.line[:150], read)
		assert a.is_proper_pair and b.is_proper_pair
		assert a.tlen == H.tlen_5prime(a, b)
		assert abs(abs(a.tlen) - p.insert) <= 10


def test_c_api_map_one(synth, tmp_path):
	"""mb_map() through a C driver: each read's SAM-primary hit is on the true contig and strand."""
	reads = H.simulate_reads(synth.genome, 60, seed=23, unique=True, qual=None)
	fa = H.write_fasta(str(tmp_path / "r.fa"), reads)
	tsv = str(tmp_path / "truth.tsv")
	with open(tsv, "w") as f:
		for r in reads:
			f.write("%s\t%s\t%d\t%d\t%s\n" % (r.name, r.contig, r.pos, r.end, r.strand))
	H.assert_c_case("drv_pilot", "map_one", synth.index, fa, tsv, 3)


@pytest.mark.external
def test_samtools_reads_the_output(chrm_pe, tmp_path):
	samtools = shutil.which("samtools")
	if samtools is None:
		pytest.skip("samtools not on PATH")
	p = str(tmp_path / "x.sam")
	with open(p, "w") as f:
		f.write(chrm_pe.stdout)
	r = subprocess.run([samtools, "view", "-c", p], stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=60)
	assert r.returncode == 0, r.stderr
	assert int(r.stdout) == len(chrm_pe.sam().records)


@pytest.mark.external
def test_pysam_reads_the_output(chrm_pe, tmp_path):
	pysam = pytest.importorskip("pysam")
	p = str(tmp_path / "x.sam")
	with open(p, "w") as f:
		f.write(chrm_pe.stdout)
	with pysam.AlignmentFile(p, "r") as f:
		recs = list(f)
	ours = chrm_pe.sam().records
	assert len(recs) == len(ours)
	assert [(r.query_name, r.flag, r.reference_start) for r in recs] == [(r.qname, r.flag, r.pos0) for r in ours]
