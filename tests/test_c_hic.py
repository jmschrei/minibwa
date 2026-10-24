"""Area C: -5 (MB_F_PRIMARY5, mb_set_sam_pri with is_primary5 and its re-selection
in mb_pair) and --hic (-5P).

minibwa's usage text: "-5 take the alignment with the smallest query position as
primary" and "--hic map Hi-C reads; equivalent to option -5P"; minibwa.1 adds
"if a read has supplementary alignments". The query position is on the read
as sequenced, so for a reverse-strand alignment it is the clip on the right
of the CIGAR. When -5 moves the primary away from the hit that was paired, no
record of the pair keeps flag 0x2.
"""

import pytest

import minibwa_helpers as H
import helpers_c as C


def _qs(rec, qlen):
	"""Query start of a record on the read as sequenced."""
	return rec.right_clip if rec.is_reverse else rec.left_clip


def _at(rec, hit, tol=3):
	"""rec aligns hit = (contig, start, reverse, length), give or take `tol` bases at either end
	(a base at a chimeric junction can match both sides by chance)."""
	c, ts, rev, n = hit
	return rec.rname == c and rec.is_reverse == rev and abs(rec.pos0 - ts) <= tol and abs(rec.ref_end0 - ts - n) <= tol


def _build():
	pg = C.PlantedGenome(701, [("chrA", 60000), ("chrB", 30000)])
	g = pg.genome()
	sc = []
	# chimeric read 1, read 2 a plain mate; (name, layout) where the part next to read 2 is
	# "pair" and the other part comes from chrB
	b = 1000
	for name, five_pairs, five_len, r1_rev in (("c5else_f", False, 50, False), ("c5pair_f", True, 100, False),
			("c5else_r", False, 50, True), ("c5pair_r", True, 100, True)):
		b += 2000
		L = 150
		pl = five_len if five_pairs else L - five_len  # bases next to read 2
		other = g["chrB"][b // 2:b // 2 + (L - pl)]
		if not r1_rev:  # pairing part forward at [b, b+pl); read 2 reverse ending at b+400
			part = g["chrA"][b:b + pl]
			r2 = C.mate_at(g, name, "chrA", b + 250, True)
			pair_hit = ("chrA", b, False, pl)
		else:  # pairing part reverse ending at b+400; read 2 forward at b
			part = H.revcomp(g["chrA"][b + 400 - pl:b + 400])
			r2 = C.mate_at(g, name, "chrA", b, False)
			pair_hit = ("chrA", b + 400 - pl, True, pl)
		seq = part + other if five_pairs else other + part
		r1 = C.Mate(name, seq, C._qual(name, L), "chrA", pair_hit[1], pair_hit[1] + pl, r1_rev, [])
		p = C.Pair(name, r1, r2)
		p.pair_hit, p.other_hit, p.five_pairs = pair_hit, ("chrB", b // 2, False, L - pl), five_pairs
		sc.append(p)
	return g, sc


@pytest.fixture(scope="module")
def chim(index_of, tmp_path_factory):
	g, sc = _build()
	idx = index_of(g)
	f1, f2 = C.write_pairs_c(str(tmp_path_factory.mktemp("chim")), sc)
	runs = dict((k, H.run_map(idx, [f1, f2], a)) for k, a in (("default", []), ("p5", ["-5"]), ("hic", ["--hic"]),
		("p5P", ["-5", "-P"]), ("P", ["-P"])))
	return g, sc, runs


def _recs(sam, name, seg):
	return [r for r in sam.records if r.qname == name and r.segment == seg]


def test_chimeric_read_has_two_primary_hits(chim):
	"""Precondition: read 1 of every scenario has one record next to read 2 and one on chrB
	(primary plus supplementary) in every mode."""
	g, sc, runs = chim
	for k, res in runs.items():
		sam = res.sam()
		assert H.pair_invariant_errors(sam) == [], k
		for p in sc:
			recs = _recs(sam, p.name, 1)
			assert len(recs) == 2 and sum(r.is_supplementary for r in recs) == 1, (k, p.name)
			assert any(_at(r, p.pair_hit) for r in recs) and any(_at(r, p.other_hit) for r in recs), (k, p.name)


@pytest.mark.parametrize("mode", ["p5", "hic", "p5P"])
def test_primary5_takes_hit_at_read_start(chim, mode):
	"""With -5 (alone, with -P, or as --hic) read 1's primary record is the hit covering its first
	base (query start 0, whatever the strand), and the other hit is supplementary."""
	g, sc, runs = chim
	sam = runs[mode].sam()
	for p in sc:
		recs = _recs(sam, p.name, 1)
		pri = [r for r in recs if r.is_primary][0]
		assert _qs(pri, 150) == 0, (mode, p.name, pri.line[:100])
		assert _at(pri, p.pair_hit if p.five_pairs else p.other_hit), (mode, p.name, pri.line[:100])
		sup = [r for r in recs if r.is_supplementary][0]
		assert _qs(sup, 150) > 0 and "H" in sup.cigar


def test_without_primary5_the_better_hit_is_primary(chim):
	"""Without -5 the primary is the higher-scoring hit, the 100 bp part next to read 2, wherever
	it lies in the read."""
	g, sc, runs = chim
	for k in ("default", "P"):
		prim = runs[k].sam().primaries()
		for p in sc:
			assert _at(prim[(p.name, 1)], p.pair_hit), (k, p.name)


def test_primary5_change_clears_proper_pair(chim):
	"""In paired mode the pair is chosen first. If -5 then moves read 1's primary off the paired
	hit, no record of the pair carries 0x2; if the paired hit is also the first in the read, 0x2
	stays on both primaries and on the supplementary record. Without -5 every primary and
	supplementary record of these pairs is proper."""
	g, sc, runs = chim
	d, p5 = runs["default"].sam(), runs["p5"].sam()
	for p in sc:
		for r in d.records:
			if r.qname == p.name:
				assert r.is_proper_pair, r.line[:80]
		recs = [r for r in p5.records if r.qname == p.name]
		assert all(r.is_proper_pair == p.five_pairs for r in recs), (p.name, [r.flag for r in recs])


def test_hic_equals_primary5_without_pairing(chim):
	"""--hic gives the same output as -5 -P, and neither sets 0x2."""
	g, sc, runs = chim
	assert H.sam_body(runs["hic"].stdout) == H.sam_body(runs["p5P"].stdout)
	assert not any(r.is_proper_pair for r in runs["hic"].sam().records)


def test_mate_fields_follow_primary5(chim):
	"""Read 2's RNEXT/PNEXT and MC:Z point at read 1's -5 primary, not at its best hit."""
	g, sc, runs = chim
	for k in ("p5", "hic"):
		prim = runs[k].sam().primaries()
		for p in sc:
			a, b = prim[(p.name, 1)], prim[(p.name, 2)]
			assert (b.rnext if b.rnext != "=" else b.rname, b.pnext) == (a.rname, a.pos), (k, p.name)
			assert b.tags["MC"] == a.cigar


# ---------------------------------------------------------------------------
# Simulated Hi-C pairs
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def hic(synth, tmp_path_factory):
	pairs = H.simulate_hic_pairs(synth.genome, 80, seed=71, chimeric_frac=0.5, junction_min=40, unique=True, unique_pad=50)
	d = tmp_path_factory.mktemp("hic")
	f1, f2 = H.write_pairs(str(d / "h_1.fq"), str(d / "h_2.fq"), pairs)
	return pairs, f1, f2, H.run_map(synth.index, [f1, f2], ["--hic"])


def _near(rec, seg, tol=5):
	return (rec.rname == seg.contig and rec.strand == seg.strand and abs(rec.pos0 - seg.pos) <= tol
		and abs(rec.ref_end0 - seg.end) <= tol)


def test_hic_mates_at_their_loci(synth, hic):
	"""Simulated Hi-C pairs (mates from independent loci, half of read 1s chimeric) with --hic:
	every primary lies at the locus of the read's first segment and covers its first base; a
	chimeric read 1 has a hard-clipped supplementary record at its second segment; no record is
	proper; SEQ, NM and the mate fields are consistent."""
	pairs, f1, f2, res = hic
	sam = res.sam()
	reads = dict((p.name, [p.r1.seq, p.r2.seq]) for p in pairs)
	assert H.sam_invariant_errors(sam, synth.genome, reads) == []
	assert H.pair_invariant_errors(sam) == []
	assert not any(r.is_proper_pair for r in sam.records)
	prim = sam.primaries()
	n_chim = 0
	for p in pairs:
		for seg, read in ((1, p.r1), (2, p.r2)):
			r = prim[(p.name, seg)]
			assert _near(r, read.segments[0]), (r.line[:100], read)
			assert _qs(r, len(read)) == 0
		supp = [r for r in sam.records if r.qname == p.name and r.is_supplementary]
		if p.r1.is_chimeric:
			n_chim += 1
			assert len(supp) == 1 and supp[0].segment == 1 and "H" in supp[0].cigar
			assert _near(supp[0], p.r1.segments[1]), (supp[0].line[:100], p.r1)
		else:
			assert supp == []
	assert n_chim >= 30


def test_hic_mates_are_mates(hic):
	"""Mates on two contigs get RNEXT set to the mate's contig and TLEN 0; mates on one contig
	the 5'-to-5' TLEN, which here exceeds 20 kbp."""
	pairs, f1, f2, res = hic
	prim = res.sam().primaries()
	n_inter = n_intra = 0
	for p in pairs:
		a, b = prim[(p.name, 1)], prim[(p.name, 2)]
		if a.rname != b.rname:
			n_inter += 1
			assert (a.rnext, a.tlen, b.rnext, b.tlen) == (b.rname, 0, a.rname, 0)
		else:
			n_intra += 1
			assert a.rnext == "=" and a.tlen == H.tlen_5prime(a, b) and abs(a.tlen) > 19000
	assert n_inter >= 20 and n_intra >= 20
