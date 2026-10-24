"""Area B: base alignments of simulated reads checked against recomputation from
the genome. Every tag and column the alignment code writes is recomputed
independently: CIGAR structure and indel placement, NM, AS (dual-affine score
of the CIGAR), ms (best-scoring segment), md, cs, ds, MD, --eqx and PAF
columns 10-11. Short reads (60-300 bp, short-read mode) and accurate long
reads (0.4-5 kbp, long-read mode) on the standard synthetic genome, with
substitutions, indels up to 8 bp and N bases."""

import pytest

import minibwa_helpers as H
import helpers_b as B


@pytest.fixture(scope="module")
def se(synth, tmp_path_factory):
	"""Three runs over one read set: SAM with cs and secondaries, SAM with MD, PAF with ds and --eqx."""
	reads = H.simulate_reads(synth.genome, 400, length=(60, 300), seed=31, sub_rate=0.02, indel_rate=0.004,
		max_indel=8, n_rate=0.003)
	fq = H.write_fastq(str(tmp_path_factory.mktemp("b_se") / "se.fq"), reads)
	return dict(reads=reads, truth=dict((r.name, r) for r in reads),
		cs=H.run_map(synth.index, fq, ["-b", "cs", "--outn", "5"]),
		md=H.run_map(synth.index, fq, ["-b", "MD"]),
		paf=H.run_map(synth.index, fq, ["-f", "-b", "ds", "--eqx", "--outn", "5"]))


@pytest.fixture(scope="module")
def lr(synth, tmp_path_factory):
	reads = H.simulate_long_reads(synth.genome, 60, length=(400, 5000), seed=41, n_rate=0.001)
	fq = H.write_fastq(str(tmp_path_factory.mktemp("b_lr") / "lr.fq"), reads)
	return dict(reads=reads, truth=dict((r.name, r) for r in reads),
		cs=H.run_map(synth.index, fq, ["-b", "cs", "--outn", "5"]),
		paf=H.run_map(synth.index, fq, ["-f", "-b", "ds", "--eqx"]))


def _records(run, genome):
	"""(record, aligned query, contig sequence) for every mapped SAM record."""
	sam = run["cs"].sam()
	out = []
	for rec in sam.records:
		if rec.is_unmapped:
			continue
		out.append((rec, B.aligned_query(rec, run["truth"][rec.qname].seq), genome[rec.rname]))
	return out


def _paf_records(run, genome):
	out = []
	for p in run["paf"].paf():
		if p.is_unmapped:
			continue
		read = run["truth"][p.qname].seq
		s = read if p.strand == "+" else H.revcomp(read)
		qs = p.qstart if p.strand == "+" else p.qlen - p.qend
		out.append((p, s[qs:qs + p.qend - p.qstart], genome[p.tname]))
	return out


@pytest.mark.parametrize("which", ["se", "lr"])
def test_sam_consistency(synth, se, lr, which):
	"""SEQ is the read, CIGAR spans it and stays on its contig, NM equals the recomputed edit
	distance; every simulated read has exactly one primary record."""
	run = se if which == "se" else lr
	sam = run["cs"].sam()
	assert H.sam_invariant_errors(sam, synth.genome, dict((r.name, [r.seq]) for r in run["reads"])) == []
	assert sum(not r.is_unmapped for r in sam.primaries().values()) >= 0.97 * len(run["reads"])


@pytest.mark.parametrize("which", ["se", "lr"])
def test_cigar_structure(synth, se, lr, which):
	"""No zero-length or repeated operations, clips only at the ends, no indel next to a
	clip (mb_fix_cigar), and every indel between matches left-aligned (it cannot move one
	base left with the same alignment)."""
	run = se if which == "se" else lr
	errs = []
	for rec, q, ref in _records(run, synth.genome):
		for e in B.cigar_problems(rec.cigar) + B.indels_left_aligned(ref, rec.pos0, q, rec.cigar):
			errs.append("%s %s: %s" % (rec.qname, rec.cigar, e))
	assert errs == []


@pytest.mark.parametrize("which", ["se", "lr"])
def test_alignment_scores(synth, se, lr, which):
	"""AS is the score of the reported alignment under the default scoring (match 2,
	mismatch -8, N -1, gap min(12+2l, 23+l)); ms is its best-scoring segment with gap cost
	q+e (short reads) or q+e*log2(1+l) (long reads); 0 <= md <= ms."""
	run = se if which == "se" else lr
	errs = []
	n = 0
	for rec, q, ref in _records(run, synth.genome):
		n += 1
		want_as = B.aln_score(ref, rec.pos0, q, rec.cigar)
		want_ms = B.segment_score(ref, rec.pos0, q, rec.cigar, long_mode=which == "lr")
		if rec.tags["AS"] != want_as:
			errs.append("%s %s: AS %d, alignment scores %d" % (rec.qname, rec.cigar, rec.tags["AS"], want_as))
		if rec.tags["ms"] != want_ms:
			errs.append("%s %s: ms %d, recomputed %d" % (rec.qname, rec.cigar, rec.tags["ms"], want_ms))
		if not 0 <= rec.tags["md"] <= rec.tags["ms"]:
			errs.append("%s: md %d outside [0, ms %d]" % (rec.qname, rec.tags["md"], rec.tags["ms"]))
	assert n > 0 and errs == []


def test_short_reads_at_least_as_good_as_truth(synth, se):
	"""For reads from unique sequence without N, the reported alignment scores at least
	as much as the simulated (true) alignment, and the read lands at its true position.
	Only short reads (< 100 bp) with several errors may stay unmapped (too few seeds)."""
	n = 0
	for rec in se["cs"].sam().primaries().values():
		t = se["truth"][rec.qname]
		if not synth.genome.is_unique(t.contig, t.pos, t.end) or "N" in t.seq:
			continue
		if rec.is_unmapped:
			assert len(t.seq) < 100 and t.nm >= 3, t
			continue
		n += 1
		truth_score = B.aln_score(synth.genome[t.contig], t.pos, t.fwd_seq, t.cigar)
		assert rec.tags["AS"] >= truth_score, (rec.line[:200], t)
		err = H.placement_error(rec, t)
		assert err is not None and err <= 8, (rec.line[:200], t)
		assert rec.tags["md"] == rec.tags["ms"]  # nothing else scores on these reads
	assert n > 200


@pytest.mark.parametrize("which", ["se", "lr"])
def test_cs(synth, se, lr, which):
	"""cs:Z equals its recomputation from the genome (bases compared as A/C/G/T/N codes)."""
	run = se if which == "se" else lr
	bad = [(rec.qname, rec.tags.get("cs"), H.compute_cs(ref, rec.pos0, q, rec.cigar)) for rec, q, ref in _records(run, synth.genome)
		if rec.tags.get("cs") != H.compute_cs(ref, rec.pos0, q, rec.cigar)]
	assert bad == []


def test_md(synth, se):
	"""MD:Z equals its recomputation (N against N counts as a match, as in cs; base minibwa
	omits the final 0 after a terminal mismatch, see known-bugs/md-trailing-zero.md)."""
	sam = se["md"].sam()
	n = 0
	for rec in sam.records:
		if rec.is_unmapped:
			continue
		n += 1
		q = B.aligned_query(rec, se["truth"][rec.qname].seq)
		md = H.compute_md(synth.genome[rec.rname], rec.pos0, q, rec.cigar, ambiguous_equal=True)
		assert rec.tags["MD"] == H.md_from_minibwa(md), (rec.line[:200], md)
		assert "cs" not in rec.tags and "ds" not in rec.tags
	assert n > 300


@pytest.mark.parametrize("which", ["se", "lr"])
def test_paf_ds_eqx_columns(synth, se, lr, which):
	"""PAF with -b ds and --eqx: cg:Z uses only =, X, I, D and equals the M-CIGAR rewritten
	from the genome (N is always X); ds:Z equals its recomputation (after merging identical
	runs that --eqx splits at N-against-N 'X' runs) and, without brackets, cs;
	columns 10-11 are matching bases and alignment length without ambiguous bases; NM is
	blen - mlen plus ambiguous bases."""
	run = se if which == "se" else lr
	n = 0
	for p, q, ref in _paf_records(run, synth.genome):
		n += 1
		cg = p.tags["cg"]
		m = B.to_m_cigar(cg)
		assert set(op for _, op in H.parse_cigar(cg)) <= set("=XID"), cg
		assert B.eqx_cigar(ref, p.tstart, q, m) == cg, (p.line[:200])
		ds = B.compute_ds(ref, p.tstart, q, m)
		assert B.merge_cs_runs(p.tags["ds"]) == ds, (p.line[:200], ds)
		assert ds.replace("[", "").replace("]", "") == H.compute_cs(ref, p.tstart, q, m)
		assert (p.nmatch, p.alen) == B.mlen_blen(ref, p.tstart, q, m), p.line[:200]
		assert p.tags["NM"] == H.compute_nm(ref, p.tstart, q, m)
	assert n > 50


def test_eqx_matches_plain_alignment(se):
	"""--eqx and -b change only the output: the alignment (positions, M-CIGAR, AS) is the same."""
	plain = {}
	for rec in se["cs"].sam().records:
		if rec.is_primary and not rec.is_unmapped:
			plain[rec.qname] = rec
	n = 0
	for p in se["paf"].paf():
		if p.is_unmapped or not p.is_primary or p.qname not in plain:
			continue
		rec = plain[p.qname]
		if rec.is_supplementary or rec.pos0 != p.tstart:
			continue
		n += 1
		core = "".join("%d%s" % (k, op) for k, op in rec.cigar_ops if op not in "SH")
		assert B.to_m_cigar(p.tags["cg"]) == core and p.tags["AS"] == rec.tags["AS"], (p.line[:150], rec.line[:150])
	assert n > 300


def test_mapq_and_secondaries(synth, se):
	"""MAPQ is 0-60, secondaries have MAPQ 0 and SEQ '*'; reads from unique sequence get MAPQ
	>= 20; reads entirely inside an exact repeat get MAPQ 0 and list the other copies."""
	sam = se["cs"].sam()
	for rec in sam.records:
		assert 0 <= rec.mapq <= 60
		if rec.is_secondary:
			assert rec.mapq == 0 and rec.seq == "*"
	n_rep = 0
	for rec in sam.primaries().values():
		t = se["truth"][rec.qname]
		if synth.genome.is_unique(t.contig, t.pos, t.end, 20) and "N" not in t.seq and not rec.is_unmapped:
			assert rec.mapq >= 20, rec.line[:200]
		feats = synth.genome.features_overlapping(t.contig, t.pos, t.end)
		if len(feats) == 1 and feats[0].kind == "repeat" and feats[0].start <= t.pos and t.end <= feats[0].end:
			n_rep += 1
			assert rec.mapq == 0, rec.line[:200]
			assert rec.tags.get("n2", 0) >= 1
	assert n_rep >= 3
