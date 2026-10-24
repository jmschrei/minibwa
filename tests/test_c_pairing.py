"""Area C: choosing the pair (mb_pair_hits, mb_pair_score) and the pair MAPQ and
paired/unpaired decision in mb_pair, with -U.

The expected outcome is recomputed in Python from the single-end alignments
that minibwa reports with -P (every hit, through --outn/--outs=0): the
candidate pairs are FR placements within the -I bounds; a pair scores the two
alignment scores (ms:i) plus .721*log(2*erfc(|z|/sqrt(2)))*A for the
normalized insert size z; the best pair is used if it scores at least the two
best single-end scores minus -U times A; the pair MAPQ follows the bwa-mem
formula that minibwa keeps. The scenarios are planted in a random genome so
that no hit has high-occurrence seeds and nothing is rescued.
"""

import math

import pytest

import minibwa_helpers as H
import helpers_c as C

AVG, STD = 600, 100
LO, HI = 200, 1000  # -I 600,100: lo = 600-400, hi = 600+400
I_OPT = "%d,%d" % (AVG, STD)
A, B, Q, E = 2, 8, 12, 2


def _build():
	pg = C.PlantedGenome(401, [("chrA", 150000), ("chrB", 50000)])
	plans = []  # (name, read 1 as (contig, 5' position, reverse), x, y): read 2 is x on the reverse strand if y, else
	# on the forward strand; or x is None and y gives read 2 as (contig, 5' position, reverse)
	b = 2000

	def block():
		nonlocal b
		b += 6000
		return b

	# A: read 2 has two or three exact copies within the bounds; the one closest to the mean wins
	for k, (anchor_rev, inserts) in enumerate([(False, (600, 800)), (False, (400, 600)), (True, (600, 800)),
			(True, (400, 600)), (False, (600, 800, 1000)), (False, (300, 600, 900)), (True, (300, 600, 900))]):
		x = pg.random(150)
		s = block()
		e = s + 1200
		for ins in inserts:
			ts = s + ins - 150 if not anchor_rev else e - ins
			pg.plant("chrA", ts, x)
		anchor = ("chrA", e if anchor_rev else s, anchor_rev)
		plans.append(("A%d" % k, anchor, x, not anchor_rev))
	# E: a fragment inside an exact 700 bp segmental duplication; two equally good pairs. A third
	# copy of read 2 on chrB has, next to it, a copy of read 1 without seeds (8 substitutions, one
	# every 19 bases): only mate rescue could find it, and rescue is skipped when the two best
	# pairs tie.
	dup = pg.random(700)
	s = block()
	pg.plant("chrA", s, dup)
	pg.plant("chrB", 3000, dup)
	pg.plant("chrB", 40000 + 600 - 150, dup[500:650])
	pg.plant("chrB", 40000, C.substitute(dup[50:200], C.unseedable_subs(step=19, first=9)))
	plans.append(("E0", ("chrA", s + 50, False), None, ("chrA", s + 50 + 600, True)))
	# R: both reads have two exact copies, but only one combination pairs; the pair MAPQ is
	# then capped by the pair score plus 4*A minus the sum of the second-best single-end scores
	x = pg.random(150)
	y = pg.random(150)
	s = block()
	pg.plant("chrA", s, y)
	pg.plant("chrB", 30000, y)
	pg.plant("chrA", s + 600 - 150, x)
	pg.plant("chrB", 35000, x)
	plans.append(("R0", ("chrA", s, False), x, True))
	# F: read 2 has one copy within the bounds and one on chrB
	x = pg.random(150)
	s = block()
	pg.plant("chrA", s + 600 - 150, x)
	pg.plant("chrB", 10000, x)
	plans.append(("F0", ("chrA", s, False), x, True))
	# C: read 2 differs from its chrB copy at one base and from the copy near read 1 at three;
	# the pair wins only if -U is large enough. In C1 the insert size (670) puts the pair score
	# exactly 2*U below the single-end scores at the threshold U.
	for k, ins in enumerate((600, 670)):
		x = pg.random(150)
		s = block()
		pg.plant("chrB", 20000 + k * 1000, C.substitute(x, [75]))
		pg.plant("chrA", s + ins - 150, C.substitute(x, [40, 80, 120]))
		plans.append(("C%d" % k, ("chrA", s, False), x, True))
	# G: a plain unique pair, and one with read 1 reverse
	s = block()
	plans.append(("G0", ("chrA", s, False), None, ("chrA", s + 500, True)))
	s = block()
	plans.append(("G1", ("chrA", s + 500, True), None, ("chrA", s, False)))
	g = pg.genome()
	pairs = []
	for name, (c, p5, rev), x, mate in plans:
		r1 = C.mate_at(g, name, c, p5 - 150 if rev else p5, rev)
		if x is not None:  # read 2 is the sequence x on the strand opposite to read 1
			seq = H.revcomp(x) if mate else x
			r2 = C.Mate(name, seq, C._qual(name + "/2", 150), None, 0, 0, mate, [])
		else:
			mc, mp5, mrev = mate
			r2 = C.mate_at(g, name, mc, mp5 - 150 if mrev else mp5, mrev)
		pairs.append(C.Pair(name, r1, r2))
	return g, pairs


@pytest.fixture(scope="module")
def pairing(index_of, tmp_path_factory):
	g, pairs = _build()
	idx = index_of(g)
	d = tmp_path_factory.mktemp("pairing")
	f1, f2 = C.write_pairs_c(str(d), pairs)
	se = H.run_map(idx, [f1, f2], ["-P", "--outn=50", "--outs=0"])
	return g, pairs, idx, f1, f2, se


class Hit(object):
	"""A single-end hit as reported by a -P run."""

	def __init__(self, rec):
		self.contig, self.ts, self.te, self.rev = rec.rname, rec.pos0, rec.ref_end0, rec.is_reverse
		self.dp = rec.tags["ms"]
		self.mapq = rec.mapq
		self.primary = not rec.is_secondary
		self.blen = sum(n for n, op in rec.cigar_ops if op in "MID")
		self.mlen = self.blen - rec.tags["NM"]

	@property
	def p5(self):
		return self.te if self.rev else self.ts

	def key(self):
		return (self.contig, self.ts, self.rev)


def _se_hits(se):
	"""{(qname, segment): [Hit]} from a -P run with every hit written."""
	out = {}
	for rec in se.sam().records:
		assert not rec.is_unmapped and not rec.is_supplementary, rec.line[:100]
		out.setdefault((rec.qname, rec.segment), []).append(Hit(rec))
	return out


def expected_pairing(h1, h2, U):
	"""(primary hit keys allowed for read 1 and read 2 as a set of pairs, proper, mapq1, mapq2)."""
	cands = []
	for x in h1:
		for y in h2:
			if x.contig != y.contig:
				continue
			left, right = (x, y) if x.p5 < y.p5 else (y, x)
			dist = right.p5 - left.p5
			if left.rev or not right.rev or not LO <= dist <= HI:
				continue
			cands.append((C.pair_score(x.dp, y.dp, dist, AVG, STD, A), x, y))
	se1 = sorted((h.dp for h in h1), reverse=True) + [0]
	se2 = sorted((h.dp for h in h2), reverse=True) + [0]
	score_se, score_se2 = se1[0] + se2[0], se1[1] + se2[1]
	pri1 = [h for h in h1 if h.primary][0]
	pri2 = [h for h in h2 if h.primary][0]
	if not cands:
		return set([(pri1.key(), pri2.key())]), False, pri1.mapq, pri2.mapq
	scores = sorted((s for s, _, _ in cands), reverse=True)
	best = scores[0]
	sub = scores[1] if len(scores) > 1 else 0
	n_sub = sum(1 for s in scores if s >= best - max(A + B, Q + E))
	if best < score_se - U * A:  # unpaired: cap every MAPQ
		cap = 6 * (score_se - U * A - best) // A
		return set([(pri1.key(), pri2.key())]), False, min(pri1.mapq, cap), min(pri2.mapq, cap)
	allowed, mapqs = set(), set()
	for s, x, y in cands:
		if s != best:
			continue
		score2 = sub
		if (not x.primary or not y.primary) and score2 < score_se - U * A:
			score2 = score_se - U * A
		diff = min(best - score2, best + 4 * A - score_se2)
		ident = float(x.mlen + y.mlen) / (x.blen + y.blen)
		mq = int(6.02 * ident * ident * diff / A - 4.343 * math.log(n_sub + 1) + .499)
		mq = max(mq, 0)
		mq = min(mq, 60)
		if mq <= 0 and best > score2:
			mq = 1
		res = []
		for h, hits in ((x, h1), (y, h2)):
			if len(hits) == 1:
				res.append(max(h.mapq, mq))
			else:
				res.append(mq if h.mapq > mq else int(.2 * h.mapq + .8 * mq + .499))
		allowed.add((x.key(), y.key()))
		mapqs.add(tuple(res))
	assert len(mapqs) == 1, "tied pairs would give different MAPQs; the scenario is ambiguous"
	m1, m2 = mapqs.pop()
	return allowed, True, m1, m2


def _check_run(pairs, se_hits, res, U):
	sam = res.sam()
	assert H.pair_invariant_errors(sam) == []
	prim = sam.primaries()
	for p in pairs:
		a, b = prim[(p.name, 1)], prim[(p.name, 2)]
		allowed, proper, m1, m2 = expected_pairing(se_hits[(p.name, 1)], se_hits[(p.name, 2)], U)
		got = ((a.rname, a.pos0, a.is_reverse), (b.rname, b.pos0, b.is_reverse))
		assert got in allowed, (U, p.name, got, allowed)
		assert a.is_proper_pair == b.is_proper_pair == proper, (U, p.name, a.flag, b.flag)
		assert (a.mapq, b.mapq) == (m1, m2), (U, p.name)
		if a.rname == b.rname:
			assert a.tlen == H.tlen_5prime(a, b)


def test_single_end_hits_as_designed(pairing):
	"""Precondition for the reference: read 2 of the A scenarios has one hit per planted copy,
	E's mates two hits each, and the G pairs one hit per mate."""
	g, pairs, idx, f1, f2, se = pairing
	hits = _se_hits(se)
	n_copies = {"A0": 2, "A1": 2, "A2": 2, "A3": 2, "A4": 3, "A5": 3, "A6": 3, "F0": 2, "C0": 2, "C1": 2}
	for name, n in n_copies.items():
		assert len(hits[(name, 2)]) == n, name
		assert len(hits[(name, 1)]) == 1, name
	assert len(hits[("E0", 1)]) == 2 and len(hits[("E0", 2)]) == 3
	assert len(hits[("R0", 1)]) == len(hits[("R0", 2)]) == 2
	for name in ("G0", "G1"):
		assert len(hits[(name, 1)]) == len(hits[(name, 2)]) == 1


@pytest.mark.parametrize("U", [17, 0, 40])
def test_pair_choice_and_mapq_match_reference(pairing, U):
	"""For every scenario, the primary of each mate, flag 0x2 and both MAPQs equal the Python
	recomputation from the single-end hits: the copy of read 2 closest to the mean insert size
	wins, tied pairs leave both mates on one copy with MAPQ 0, and an unpaired outcome caps the
	MAPQs at 6*(score deficit)/A."""
	g, pairs, idx, f1, f2, se = pairing
	res = H.run_map(idx, [f1, f2], ["-I", I_OPT, "-U", U])
	_check_run(pairs, _se_hits(se), res, U)


def test_closest_copy_to_mean_insert_wins(pairing):
	"""Independently of the reference: read 2's primary is the planted copy whose insert size is
	closest to the -I mean, whichever copy single-end mapping made primary."""
	g, pairs, idx, f1, f2, se = pairing
	res = H.run_map(idx, [f1, f2], ["-I", I_OPT])
	prim = res.sam().primaries()
	for p in pairs:
		if not p.name.startswith("A"):
			continue
		a, b = prim[(p.name, 1)], prim[(p.name, 2)]
		ins = abs((b.ref_end0 if b.is_reverse else b.pos0) - (a.ref_end0 if a.is_reverse else a.pos0))
		assert ins == AVG and a.is_proper_pair, (p.name, a.line[:80], b.line[:80])


def test_tied_pairs_keep_mates_together(pairing):
	"""A fragment from an exact segmental duplication pairs on one copy (both mates on the
	same contig, proper, MAPQ 0). With two equally good pairs minibwa skips mate rescue (one
	of the README's "heuristics to skip unnecessary mate rescue"): read 1 gains no hit next
	to read 2's third, unpaired copy, where only rescue could place it. -U 50 lowers the
	rescue's minimum score (sum of read 1's scores / A - U = 100) below that copy's local score
	(142 matches - 8 * 4 = 110), so it is the skip, not the score, that keeps it out."""
	g, pairs, idx, f1, f2, se = pairing
	res = H.run_map(idx, [f1, f2], ["-I", I_OPT, "-U", "50", "--outn=50", "--outs=0"])
	prim = res.sam().primaries()
	a, b = prim[("E0", 1)], prim[("E0", 2)]
	assert a.rname == b.rname and a.is_proper_pair and b.is_proper_pair
	assert (a.mapq, b.mapq) == (0, 0)
	assert abs(a.tlen) == 600
	r1 = [r for r in res.sam().records if r.qname == "E0" and r.segment == 1]
	assert sorted((r.rname, r.pos0) for r in r1) == sorted((h.contig, h.ts) for h in _se_hits(se)[("E0", 1)])


def test_unpair_penalty_threshold(pairing):
	"""-U sets the penalty for leaving a pair unpaired. Scenario C pairs only when the pair score
	is at least the best single-end scores minus U*A: below the threshold U* (from the
	single-end scores) read 2 stays on its better chrB hit without 0x2, at U* it moves next to
	read 1 with 0x2. In C1 the pair score equals the single-end scores minus U*A exactly."""
	g, pairs, idx, f1, f2, se = pairing
	hits = _se_hits(se)
	n_equal = 0
	for name in ("C0", "C1"):
		h1, h2 = hits[(name, 1)], hits[(name, 2)]
		near = [h for h in h2 if h.contig == "chrA"][0]
		far = [h for h in h2 if h.contig == "chrB"][0]
		best = C.pair_score(h1[0].dp, near.dp, near.te - h1[0].ts, AVG, STD, A)
		score_se = h1[0].dp + far.dp
		u_star = -(-(score_se - best) // A)
		assert 2 <= u_star <= 30
		n_equal += best == score_se - u_star * A
		pc = [p for p in pairs if p.name == name]
		for U, paired in ((u_star - 1, False), (u_star, True)):
			res = H.run_map(idx, [f1, f2], ["-I", I_OPT, "-U", U])
			_check_run(pc, hits, res, U)
			prim = res.sam().primaries()
			b = prim[(name, 2)]
			assert (b.rname == "chrA") == paired and b.is_proper_pair == paired, (name, U, b.line[:100])
	assert n_equal == 1
