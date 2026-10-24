"""Area C: mate rescue (mb_matesw, mb_matesw_core, mb_matesw_align, mb_ungap,
mb_hit_sum_score) and the options that control it: --rescue, -U, -P, --pe, -I,
--dbg-aln-pe.

A mate is made "unseedable" by a substitution every 15 bases, so it has no
exact match of 19 bases (-k) and no single-end hit: any alignment it gets
comes from rescue next to its seeded mate. Rescue looks for the mate where the
insert-size bounds allow it: for a forward anchor whose 5' end is at s, a
reverse mate within [s+lo-len, s+hi); for a reverse anchor ending at e, a
forward mate within [e-hi, e-lo+len).
"""

import math

import pytest

import minibwa_helpers as H
import helpers_c as C

LO, HI = 300, 700  # -I 500,50
I_OPT = ["-I", "500,50"]
UNSEED = C.unseedable_subs()


def _window(anchor):
	"""[start, end) of the rescue window for the mate of a seeded anchor Mate."""
	if not anchor.rev:
		return anchor.p5 + LO - C.RLEN, anchor.p5 + HI
	return anchor.p5 - HI, anchor.p5 - LO + C.RLEN


def _build():
	pg = C.PlantedGenome(501, [("chrA", 200000), ("chrB", 60000)])
	scen = {}
	b = [1000]

	def block(n=2000):
		b[0] += n
		return b[0]

	# 1. window geometry: anchor read 1 or 2, forward or reverse, mate unseedable
	window = []
	for anchor_is_r1 in (True, False):
		for anchor_rev in (False, True):
			for ins in (LO, HI, HI + 1, LO - 1, HI + 200):
				window.append((anchor_is_r1, anchor_rev, ins, block()))
	scen["window"] = window
	# 2. --rescue=N: read 1 has three exact copies; the unseedable mate lies next to copy 1 only
	y = pg.random(150)
	copies = [block(), block(), block()]
	for s in copies:
		pg.plant("chrA", s, y)
	scen["ncand"] = (y, copies)
	# 3. -U and rescue candidates: read 1 differs from its chrB copy at one base and from a weak
	# chrA copy at four; the unseedable mate lies next to the weak copy
	z = pg.random(150)
	pg.plant("chrB", 5000, C.substitute(z, [75]))
	s = block()
	pg.plant("chrA", s, C.substitute(z, [40, 80, 110, 130]))
	scen["ucand"] = (z, s)
	# 4. rescue score threshold: read 2 maps perfectly elsewhere (one hit, or as a chimera of two
	# 75 bp parts); a copy with 4 substitutions lies next to read 1
	x = pg.random(150)
	pg.plant("chrB", 15000, x)
	s1 = block()
	pg.plant("chrA", s1 + 400 - 150, C.substitute(x, [30, 60, 90, 120]))
	xc = pg.random(150)
	pg.plant("chrB", 25000, xc[:75])
	pg.plant("chrB", 35000, xc[75:])
	s2 = block()
	pg.plant("chrA", s2 + 400 - 150, C.substitute(xc, [30, 60, 90, 120]))
	scen["minsc"] = ((x, s1), (xc, s2))
	# 5. 7-mer filter: a mate with a substitution every 7 bases has no shared 7-mer, one with a
	# substitution every 9 bases has many
	scen["hough"] = (block(), block())
	return pg.genome(), scen


def _window_pairs(g, scen):
	pairs = []
	for anchor_is_r1, anchor_rev, ins, s in scen["window"]:
		name = "w%s%s%d" % ("1" if anchor_is_r1 else "2", "r" if anchor_rev else "f", ins)
		if not anchor_rev:
			anchor = C.mate_at(g, name, "chrA", s, False)
			mate = C.mate_at(g, name, "chrA", s + ins - 150, True, subs=UNSEED)
		else:
			e = s + 1200
			anchor = C.mate_at(g, name, "chrA", e - 150, True)
			mate = C.mate_at(g, name, "chrA", e - ins, False, subs=UNSEED)
		p = C.Pair(name, anchor, mate) if anchor_is_r1 else C.Pair(name, mate, anchor)
		p.anchor, p.mate, p.ins = anchor, mate, ins
		pairs.append(p)
	return pairs


@pytest.fixture(scope="module")
def rescue_genome(index_of):
	g, scen = _build()
	return g, scen, index_of(g)


@pytest.fixture(scope="module")
def window_run(rescue_genome, tmp_path_factory):
	g, scen, idx = rescue_genome
	pairs = _window_pairs(g, scen)
	f1, f2 = C.write_pairs_c(str(tmp_path_factory.mktemp("window")), pairs)
	return pairs, f1, f2, H.run_map(idx, [f1, f2], I_OPT)


def _mate_recs(sam, p):
	seg = 2 if p.mate is p.r2 else 1
	recs = [r for r in sam.records if r.qname == p.name and r.segment == seg]
	arecs = [r for r in sam.records if r.qname == p.name and r.segment == 3 - seg]
	return recs, arecs


def test_unseedable_mates_have_no_single_end_hit(rescue_genome, window_run):
	"""Precondition: with pairing off (-P) every unseedable mate is unmapped and every anchor
	maps at its truth."""
	g, scen, idx = rescue_genome
	pairs, f1, f2, res = window_run
	sam = H.run_map(idx, [f1, f2], ["-P"]).sam()
	for p in pairs:
		recs, arecs = _mate_recs(sam, p)
		assert len(recs) == 1 and recs[0].is_unmapped
		a = [r for r in arecs if r.is_primary][0]
		assert (a.rname, a.pos0, a.strand, a.cigar) == ("chrA", p.anchor.ts, p.anchor.strand, "150M")


def test_mate_rescued_at_both_insert_bounds(rescue_genome, window_run):
	"""A mate whose insert size is exactly lo or hi is rescued in full: at its true position and
	strand, CIGAR 150M, NM equal to its substitutions, both mates proper, TLEN +-insert. This
	holds whichever read is the anchor and on either strand."""
	g, scen, idx = rescue_genome
	pairs, f1, f2, res = window_run
	sam = res.sam()
	assert H.pair_invariant_errors(sam) == []
	n = 0
	for p in pairs:
		if p.ins not in (LO, HI):
			continue
		recs, arecs = _mate_recs(sam, p)
		assert len(recs) == 1, p
		m = recs[0]
		assert (m.rname, m.pos0, m.strand, m.cigar) == ("chrA", p.mate.ts, p.mate.strand, "150M"), (p, m.line[:120])
		assert m.tags["NM"] == H.compute_nm(g["chrA"], m.pos0, m.seq, m.cigar) == 10
		a = [r for r in arecs if r.is_primary][0]
		assert m.is_proper_pair and a.is_proper_pair
		assert abs(m.tlen) == p.ins and m.tlen == H.tlen_5prime(m, a)
		n += 1
	assert n == 8


def test_rescued_alignment_stays_inside_window(rescue_genome, window_run):
	"""A mate one base beyond the bounds may be rescued only partly: any alignment it gets lies
	inside the window. A mate entirely outside the window is not rescued: it is unmapped and
	placed at its mate."""
	g, scen, idx = rescue_genome
	pairs, f1, f2, res = window_run
	sam = res.sam()
	n_mapped = 0
	for p in pairs:
		recs, arecs = _mate_recs(sam, p)
		a = [r for r in arecs if r.is_primary][0]
		assert (a.rname, a.pos0, a.cigar) == ("chrA", p.anchor.ts, "150M")
		m = recs[0]
		if p.ins == HI + 200:
			assert len(recs) == 1 and m.is_unmapped and (m.rname, m.pos) == (a.rname, a.pos)
			assert not a.is_proper_pair
			continue
		if m.is_unmapped:
			continue
		w0, w1 = _window(p.anchor)
		assert w0 <= m.pos0 and m.ref_end0 <= w1, (p, m.line[:120], (w0, w1))
		n_mapped += 1
	assert n_mapped >= 9  # 8 at the bounds and at least one partial rescue


@pytest.mark.parametrize("args", [["-P"], ["--rescue=0"], ["--pe=no"]], ids=["P", "rescue0", "pe_no"])
def test_rescue_disabled(rescue_genome, window_run, args):
	"""-P, --rescue=0 and --pe=no (two read files) rescue nothing: every unseedable mate stays
	unmapped and no pair is proper. As no mate has a hit of its own, all three give the same
	output."""
	g, scen, idx = rescue_genome
	pairs, f1, f2, res = window_run
	out = H.run_map(idx, [f1, f2], I_OPT + args)
	sam = out.sam()
	assert H.pair_invariant_errors(sam) == []
	assert not any(r.is_proper_pair for r in sam.records)
	for p in pairs:
		recs, _ = _mate_recs(sam, p)
		assert [r.is_unmapped for r in recs] == [True]
	if args != ["-P"]:  # no mate has a hit of its own, so nothing differs from -P
		ref = H.run_map(idx, [f1, f2], I_OPT + ["-P"])
		assert H.sam_body(out.stdout) == H.sam_body(ref.stdout)


def test_rescue_candidates_limit(rescue_genome, tmp_path):
	"""--rescue=N tries at most N single-end hits of the anchor, in hit order. Read 1 has three
	exact copies and its unseedable mate lies next to one of them; twelve pairs with the same
	sequences and different names put that copy at different ranks (names break ties). The mate
	is rescued iff the copy is among the first N hits, as listed by a -P --outn run; once
	rescued, read 1 moves to that copy and the pair is proper."""
	g, scen, idx = rescue_genome
	y, copies = scen["ncand"]
	pairs = []
	for i in range(12):
		name = "n%d" % i
		r1 = C.Mate(name, y, C._qual(name, 150), "chrA", copies[1], copies[1] + 150, False, [])
		r2 = C.mate_at(g, name, "chrA", copies[1] + 400 - 150, True, subs=UNSEED)
		pairs.append(C.Pair(name, r1, r2))
	f1, f2 = C.write_pairs_c(str(tmp_path), pairs)
	se = H.run_map(idx, [f1, f2], ["-P", "--outn=50", "--outs=0"]).sam()
	rank = {}
	for p in pairs:
		recs = [r for r in se.records if r.qname == p.name and r.segment == 1]
		assert sorted(r.pos0 for r in recs) == copies
		rank[p.name] = [r.pos0 for r in recs].index(copies[1])
	assert len(set(rank.values())) == 3  # precondition: the true copy appears at every rank
	for n in (0, 1, 2, 3):
		sam = H.run_map(idx, [f1, f2], ["--rescue=%d" % n]).sam()
		prim = sam.primaries()
		for p in pairs:
			a, b = prim[(p.name, 1)], prim[(p.name, 2)]
			want = rank[p.name] < n
			assert (not b.is_unmapped) == want, (n, p.name, rank[p.name])
			if want:
				assert (a.pos0, b.pos0, b.strand) == (copies[1], p.r2.ts, "-") and a.is_proper_pair and b.is_proper_pair


def test_unpair_penalty_limits_rescue_candidates(rescue_genome, tmp_path):
	"""Only hits scoring at least the anchor's best minus U*A are rescue candidates. Read 1's best
	hit is on chrB; its weak chrA copy, next to the unseedable mate, becomes a candidate at
	U* = ceil((best - weak)/A) (from the single-end scores), not at U*-1. When the mate is
	rescued the pair wins and read 1 moves to the weak copy."""
	g, scen, idx = rescue_genome
	z, s = scen["ucand"]
	name = "uc"
	r1 = C.Mate(name, z, C._qual(name, 150), "chrA", s, s + 150, False, [])
	r2 = C.mate_at(g, name, "chrA", s + 400 - 150, True, subs=UNSEED)
	f1, f2 = C.write_pairs_c(str(tmp_path), [C.Pair(name, r1, r2)])
	se = H.run_map(idx, [f1, f2], ["-P", "--outn=50", "--outs=0"]).sam()
	hits = dict((r.rname, r.tags["ms"]) for r in se.records if r.segment == 1)
	assert set(hits) == set(["chrA", "chrB"]) and hits["chrB"] > hits["chrA"]
	u_star = -(-(hits["chrB"] - hits["chrA"]) // 2)
	for U, rescued in ((u_star - 1, False), (u_star, True)):
		prim = H.run_map(idx, [f1, f2], ["-U", U]).sam().primaries()
		a, b = prim[(name, 1)], prim[(name, 2)]
		assert (not b.is_unmapped) == rescued, U
		if rescued:
			assert (a.rname, a.pos0, b.pos0) == ("chrA", s, s + 250) and a.is_proper_pair
		else:
			assert a.rname == "chrB" and not a.is_proper_pair


def _primary_spans(recs, qlen):
	"""(qs, qe, dp) on the read's forward strand for each primary or supplementary record."""
	out = []
	for r in recs:
		if r.is_secondary or r.is_unmapped:
			continue
		lc, rc = r.left_clip, r.right_clip
		qs, qe = (rc, qlen - lc) if r.is_reverse else (lc, qlen - rc)
		out.append((qs, qe, r.tags["ms"]))
	return out


def _sum_score(spans):
	"""Score of a read's primary hits with query overlaps counted once (overlapping parts
	prorated), as the rescue's minimum score uses it."""
	sc = qe = 0
	for s, e, dp in sorted(spans):
		if e <= qe:
			continue
		if s < qe:
			sc += int(float(e - qe) / (e - s) * dp + .499)
		else:
			sc += dp
		qe = e
	return sc


@pytest.mark.parametrize("which", [0, 1], ids=["one_hit", "chimera"])
def test_rescue_needs_score_near_own_hits(rescue_genome, tmp_path, which):
	"""A rescued alignment must score (in units of -A, mismatch ceil(-B/-A)) at least the sum of
	the mate's own primary hits, divided by -A, minus -U. Read 2 maps perfectly elsewhere, as one
	hit or as a chimera of two 75 bp parts; its copy next to read 1 has 4 substitutions and
	scores 146 - 4*4 = 130. Rescue happens from U = sum/A - 130 on, not below; the chimera's
	parts add up like one hit."""
	g, scen, idx = rescue_genome
	x, s = scen["minsc"][which]
	name = "ms%d" % which
	r1 = C.mate_at(g, name, "chrA", s, False)
	r2 = C.Mate(name, H.revcomp(x), C._qual(name, 150), None, 0, 0, True, [])
	f1, f2 = C.write_pairs_c(str(tmp_path), [C.Pair(name, r1, r2)])
	se = H.run_map(idx, [f1, f2], ["-P"]).sam()
	recs = [r for r in se.records if r.segment == 2]
	assert all(r.rname == "chrB" for r in recs) and len(recs) == which + 1
	total = _sum_score(_primary_spans(recs, 150))
	assert 290 <= total <= 310
	u_star = total // 2 - 130
	for U, rescued in ((u_star - 1, False), (u_star, True)):
		recs = [r for r in H.run_map(idx, [f1, f2], ["-U", U]).sam().records if r.segment == 2]
		near = [r for r in recs if r.rname == "chrA"]
		assert bool(near) == rescued, (U, [r.line[:80] for r in recs])
		if rescued:
			assert near[0].is_primary and near[0].is_proper_pair and near[0].pos0 == s + 250


def test_kmer_filter_before_rescue(rescue_genome, tmp_path):
	"""Rescue first looks for shared 7-mers on one diagonal. A mate with a substitution every 7
	bases shares none and is not rescued, though a local alignment would score above the
	minimum; one with a substitution every 9 bases is rescued at its true locus."""
	g, scen, idx = rescue_genome
	s7, s9 = scen["hough"]
	p7 = C.fr_pair(g, "h7", "chrA", s7, 400, subs2=list(range(3, 150, 7)))
	p9 = C.fr_pair(g, "h9", "chrA", s9, 400, subs2=list(range(4, 150, 9)))
	f1, f2 = C.write_pairs_c(str(tmp_path), [p7, p9])
	prim = H.run_map(idx, [f1, f2]).sam().primaries()
	assert prim[("h7", 2)].is_unmapped
	m = prim[("h9", 2)]
	w0, w1 = s9 + 50 - 150, s9 + 800  # the window of a forward anchor under the default bounds (50, 800)
	assert not m.is_unmapped and m.rname == "chrA" and m.is_reverse and w0 <= m.pos0 and m.ref_end0 <= w1
	assert abs(m.unclipped_end0 - p9.r2.te) <= 3


def test_dbg_aln_pe_reports_rescue(rescue_genome, window_run, tmp_path):
	"""--dbg-aln-pe (hidden) traces mate rescue on stderr and leaves stdout unchanged. The trace
	is an undocumented debugging aid, so its wording, line count and the number of rescue attempts
	it lists are not asserted: only that a pair whose mate is rescued adds trace lines to stderr."""
	g, scen, idx = rescue_genome
	pairs, f1, f2, res = window_run
	p = [q for q in pairs if q.name == "w1f%d" % LO][0]
	g1, g2 = C.write_pairs_c(str(tmp_path), [p])
	plain = H.run_map(idx, [g1, g2], I_OPT)
	dbg = H.run_map(idx, [g1, g2], I_OPT + ["--dbg-aln-pe"])
	assert H.sam_body(dbg.stdout) == H.sam_body(plain.stdout)
	m = [r for r in plain.sam().records if r.segment == 2][0]
	assert not m.is_unmapped and m.pos0 == p.mate.ts  # precondition: the mate is rescued
	trace = lambda e: [l for l in e.splitlines() if not l.startswith("[M::")]
	assert len(trace(dbg.stderr)) > len(trace(plain.stderr))


def _high_cov_build():
	"""Read 1 comes from an exact 2-copy duplication (D1, D2). Read 2 has no seed at its true
	locus T next to D1 (a substitution every 15 bases); a copy T' with one more substitution
	lies next to D2. Read 2's last 60 bases (as sequenced) occur exactly at five places on chrB,
	more than -c 3, so they count as high-occurrence sequence."""
	pg = C.PlantedGenome(502, [("chrA", 40000), ("chrB", 20000)])
	d = pg.random(150)
	t = pg.random(150)
	r = C.substitute(t, UNSEED)  # read 2 on the forward strand
	pg.plant("chrA", 5000, d)
	pg.plant("chrA", 20000, d)
	pg.plant("chrA", 5250, t)
	pg.plant("chrA", 20250, C.substitute(t, [100]))
	other = "A" if r[60] != "A" else "C"  # the base after each copy breaks the match
	for k in range(5):
		pg.plant("chrB", 2000 + 3000 * k, r[:60] + other)
	g = pg.genome()
	r1 = C.Mate("hc", d, C._qual("hc1", 150), "chrA", 5000, 5150, False, [])
	r2 = C.Mate("hc", H.revcomp(r), C._qual("hc2", 150), "chrA", 5250, 5400, True, UNSEED)
	return g, C.Pair("hc", r1, r2)


def test_rescued_hit_shares_high_occurrence_fraction(index_of, tmp_path):
	"""The pair MAPQ is scaled by 1 - (f1 + f2)/2, with f the fraction of each read covered by
	seeds occurring more than -c times. A rescued hit counts the same fraction as its read's
	other hits (60 of 150 bases here) although rescue found it without seeds. Read 1's MAPQ
	equals the recomputation with f2 = int(255*60/150)/255, from the pair scores of (D1, T) and
	(D2, T') and the alignments in the output."""
	g, p = _high_cov_build()
	idx = index_of(g)
	f1, f2 = C.write_pairs_c(str(tmp_path), [p])
	res = H.run_map(idx, [f1, f2], ["-c", "3", "-I", "400,40", "--outn=50", "--outs=0"])
	sam = res.sam()
	recs1 = [r for r in sam.records if r.segment == 1]
	recs2 = [r for r in sam.records if r.segment == 2]
	a = [r for r in recs1 if r.is_primary][0]
	b = [r for r in recs2 if r.is_primary][0]
	assert (a.rname, a.pos0, b.rname, b.pos0, b.cigar) == ("chrA", 5000, "chrA", 5250, "150M")
	assert a.is_proper_pair and b.is_proper_pair and b.tags["NM"] == 10
	d2 = [r for r in recs1 if r.pos0 == 20000]
	t2 = [r for r in recs2 if r.rname == "chrA" and r.pos0 == 20250]
	assert len(d2) == 1 and len(t2) == 1 and sum(r.rname == "chrB" for r in recs2) >= 1
	dp = lambda r: r.tags["ms"]
	s1 = C.pair_score(dp(a), dp(b), b.ref_end0 - a.pos0, 400, 40)
	s2 = C.pair_score(dp(d2[0]), dp(t2[0]), t2[0].ref_end0 - d2[0].pos0, 400, 40)
	assert s1 > s2 >= s1 - 14
	diff = min(s1 - s2, s1 + 8 - (dp(d2[0]) + dp(t2[0])))
	blen = lambda r: sum(n for n, op in r.cigar_ops if op in "MID")
	ident = float(2 * 150 - a.tags["NM"] - b.tags["NM"]) / (blen(a) + blen(b))
	mq = max(0, int(6.02 * ident * ident * diff / 2 - 4.343 * math.log(3) + .499))
	f2 = int(255. * 60 / 150)
	mq = min(60, int(mq * (1. - .5 * (0 + f2 / 255.)) + .499))
	assert 5 <= mq < 60
	assert a.mapq == int(.8 * mq + .499)


def test_rescue_window_with_lowest_bound(rescue_genome, tmp_path):
	"""-I 100,50 gives lo = max(100-200, 1) = 1, so the window of a forward anchor at s starts at
	s+1-len = s-149. An unseedable mate at [s-200, s-50) (facing away from the anchor) starts before
	the window and stays unmapped: rescue's 7-mer filter only counts diagonals on which the mate
	starts inside the window. Had lo not been raised to 1 (-100), the window would start at s-250 and
	hold the whole mate. A mate at insert 1 (touching the anchor's 5' base) is rescued."""
	g, scen, idx = rescue_genome
	s7, s9 = scen["hough"]
	s = s9 + 1000
	out = C.pair_at(g, "lo_out", "chrA", s, False, "chrA", s - 50, True, subs2=UNSEED)
	one = C.fr_pair(g, "lo_one", "chrA", s + 600, 1, subs2=UNSEED)
	f1, f2 = C.write_pairs_c(str(tmp_path), [out, one])
	prim = H.run_map(idx, [f1, f2], ["-I", "100,50"]).sam().primaries()
	assert prim[("lo_out", 2)].is_unmapped
	m = prim[("lo_one", 2)]
	assert (m.rname, m.pos0, m.strand, m.cigar) == ("chrA", one.r2.ts, "-", "150M")
