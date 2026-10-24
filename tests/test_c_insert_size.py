"""Area C: insert-size statistics (mb_pestat, mb_insert_dir, mb_select_unique_se),
the switch between inferred and fixed statistics in worker_pipeline step 1, and
the -I option (set_ins_size).

Expected statistics are recomputed in Python from the exact geometry of the
pairs (helpers_c.pestat follows the bwa-mem method minibwa documents as
"inferred"); expected proper-pair flags follow from those bounds. Tests that
check the bounds through flag 0x2 run with --rescue=0, because mate rescue can
re-align a mate that lies just outside the bounds with a clip and so make the
pair proper.
"""

import os

import pytest

import minibwa_helpers as H
import helpers_c as C


def _write(tmp_path, pairs, tag="r"):
	return C.write_pairs_c(str(tmp_path), pairs, tag)


def _flags_by_pair(sam, pairs):
	prim = sam.primaries()
	return dict((p.name, (prim[(p.name, 1)], prim[(p.name, 2)])) for p in pairs)


# ---------------------------------------------------------------------------
# Inferred statistics in four orientations, with excluded pairs
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def ins_genome(index_of):
	"""Two random contigs and an exact 300 bp repeat with two copies 6 kbp apart on chrA."""
	pg = C.PlantedGenome(301, [("chrA", 150000), ("chrB", 60000)])
	rep = pg.random(300)
	pg.plant("chrA", 140000, rep)
	pg.plant("chrA", 146000, rep)
	g = pg.genome()
	return g, index_of(g)


def _mixed_pairs(g):
	"""Unique pairs in FR (60), RF (25) and FF (3) orientation, plus pairs that must not count
	as candidates: a mate in the exact repeat (MAPQ 0; either copy would be in range), mates on
	two contigs (once at nearby coordinates), a chimeric read 1 (two primary hits), a distance
	of exactly max_pe_ins (10000), and two random mates. One FR pair at distance 9999 counts
	(outlier)."""
	rng = H.Rng(77)
	counted, other = [], []
	for i in range(60):
		ins = rng.normal_int(350, 30)
		counted.append(C.fr_pair(g, "fr%d" % i, "chrA", 1000 + i * 1500, ins, r1_rev=bool(i & 1)))
	for i in range(25):
		d = rng.normal_int(250, 20)
		left = 92000 + i * 1500
		if i & 1:  # read 1 reverse on the left, read 2 forward on the right
			counted.append(C.pair_at(g, "rf%d" % i, "chrA", left + 150, True, "chrA", left + 150 + d, False))
		else:  # read 1 forward on the right, read 2 reverse on the left
			counted.append(C.pair_at(g, "rf%d" % i, "chrA", left + 150 + d, False, "chrA", left + 150, True))
	for i in range(3):
		left = 1000 + i * 3000
		counted.append(C.pair_at(g, "ff%d" % i, "chrB", left, False, "chrB", left + 400, False))
	counted.append(C.fr_pair(g, "far9999", "chrB", 12000, 9999))
	other.append(C.fr_pair(g, "far10000", "chrB", 24000, 10000))
	# read 2 inside the repeat: two equally good hits, SE MAPQ 0
	other.append(C.pair_at(g, "rep", "chrA", 139700, False, "chrA", 140000 + 250, True))
	other.append(C.pair_at(g, "inter", "chrA", 60500, False, "chrB", 40000, True))
	other.append(C.pair_at(g, "inter_near", "chrA", 30500, False, "chrB", 30900, True))
	# chimeric read 1: 75 bp from chrA, 75 bp from chrB
	p = C.fr_pair(g, "chim", "chrA", 75000, 400)
	seq = g["chrA"][75000:75075] + g["chrB"][36000:36075]
	p.r1 = C.Mate("chim", seq, p.r1.qual, "chrA", 75000, 75075, False, [])
	other.append(p)
	junk = C.PlantedGenome(999, [("x", 300)])
	other.append(C.Pair("junk", C.Mate("junk", junk["x"][:150], "I" * 150, None, 0, 0, False, []),
		C.Mate("junk", junk["x"][150:], "I" * 150, None, 0, 0, False, [])))
	return counted, other


@pytest.fixture(scope="module")
def mixed_run(ins_genome, tmp_path_factory):
	g, idx = ins_genome
	counted, other = _mixed_pairs(g)
	pairs = counted + other
	d = tmp_path_factory.mktemp("mixed")
	f1, f2 = _write(d, pairs)
	return g, counted, other, H.run_map(idx, [f1, f2])


def _expected_stats(counted):
	dists = [[], [], [], []]
	for p in counted:
		code, dist = p.dir_dist
		dists[code].append(dist)
	return dists, C.pestat(dists)


def test_candidate_counts_per_orientation(mixed_run):
	"""The stderr report counts, per orientation, the pairs whose mates each have one primary
	hit with MAPQ >= 10, on one contig, less than 10000 bp apart (5'-to-5'). Pairs with a
	repeat mate, mates on two contigs, a chimeric mate, a distance of 10000 or an unmapped
	mate are not counted."""
	g, counted, other, res = mixed_run
	dists, _ = _expected_stats(counted)
	cands, _, _ = C.stderr_pestat(res.stderr)
	assert cands == [tuple(len(x) for x in dists)] == [(3, 61, 25, 0)]


def test_inferred_statistics_match_reference(mixed_run):
	"""Quartiles, mean, standard deviation, proper-pair bounds and the share of pairs beyond
	them equal a Python recomputation, for each orientation with at least 20 candidates and at
	least 5% of the largest orientation's count; failed orientations print nothing."""
	g, counted, other, res = mixed_run
	_, ps = _expected_stats(counted)
	_, pct, bnd = C.stderr_pestat(res.stderr)
	want_pct, want_bnd = [], []
	for d in range(4):
		s = ps[d]
		if s["failed"]:
			continue
		want_pct.append((C.DIR_NAMES[d], (s["p25"], s["p50"], s["p75"], "%.2f" % s["avg"], "%.2f" % s["std"])))
		want_bnd.append((C.DIR_NAMES[d], (s["lo"], s["hi"], "%.2f" % s["pct_lo"], "%.2f" % s["pct_hi"])))
	assert [x[0] for x in want_pct] == ["FR", "RF"]
	assert pct == want_pct
	assert bnd == want_bnd
	assert ps[1]["pct_hi"] > 0  # the 9999 bp pair is counted but lies beyond the bounds


def test_proper_flag_follows_inferred_bounds(mixed_run):
	"""Each counted pair is proper (0x2 on both primaries) iff its orientation did not fail and
	its distance lies within [lo, hi]; TLEN is the 5'-to-5' length. FR and RF pairs are proper,
	FF pairs (3 candidates) and the 9999 bp pair are not."""
	g, counted, other, res = mixed_run
	_, ps = _expected_stats(counted)
	sam = res.sam()
	assert H.pair_invariant_errors(sam) == []
	by = _flags_by_pair(sam, counted)
	n_proper = 0
	for p in counted:
		a, b = by[p.name]
		for rec, m in ((a, p.r1), (b, p.r2)):
			assert (rec.rname, rec.pos0, rec.strand, rec.cigar) == (m.contig, m.ts, m.strand, "150M"), (p, rec.line[:120])
		code, dist = p.dir_dist
		s = ps[code]
		want = not s["failed"] and s["lo"] <= dist <= s["hi"]
		assert a.is_proper_pair == b.is_proper_pair == want, (p, a.flag, b.flag)
		assert a.tlen == H.tlen_5prime(a, b) and b.tlen == -a.tlen
		n_proper += want
	assert n_proper == 85


# ---------------------------------------------------------------------------
# When statistics are inferred: 20 pairs in a batch, 20 unique pairs per orientation
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def long_inserts(ins_genome):
	"""20 FR pairs with fragments of 1100 +- 5 bp, beyond the default upper bound (800) and
	far enough beyond it that a rescue window cannot reach the mate."""
	g, idx = ins_genome
	rng = H.Rng(78)
	return [C.fr_pair(g, "L%d" % i, "chrA", 2000 + i * 2500, rng.normal_int(1100, 3), r1_rev=bool(i & 1)) for i in range(20)]


def test_fewer_than_20_pairs_use_default_statistics(ins_genome, long_inserts, tmp_path):
	"""A batch with 19 pairs prints no statistics and uses the defaults (400, 100, 800, 50):
	1100 bp fragments are not proper, and the output equals that of -I 400,100,800,50."""
	g, idx = ins_genome
	f1, f2 = _write(tmp_path, long_inserts[:19])
	res = H.run_map(idx, [f1, f2])
	assert not C.has_pestat_report(res.stderr)
	sam = res.sam()
	assert not any(r.is_proper_pair for r in sam.records)
	assert all(not r.is_unmapped for r in sam.records)
	fixed = H.run_map(idx, [f1, f2], ["-I", "400,100,800,50"])
	assert H.sam_body(res.stdout) == H.sam_body(fixed.stdout)


def test_20_pairs_infer_statistics(ins_genome, long_inserts, tmp_path):
	"""With 20 pairs the batch infers its own statistics, and the 1100 bp pairs become proper."""
	g, idx = ins_genome
	f1, f2 = _write(tmp_path, long_inserts)
	res = H.run_map(idx, [f1, f2])
	cands, pct, bnd = C.stderr_pestat(res.stderr)
	assert cands == [(0, 20, 0, 0)]
	_, ps = _expected_stats(long_inserts)
	assert bnd == [("FR", (ps[1]["lo"], ps[1]["hi"], "0.00", "0.00"))]
	sam = res.sam()
	assert sum(r.is_proper_pair for r in sam.records) == 40


def test_20_pairs_with_19_unique_disable_pairing(ins_genome, long_inserts, tmp_path):
	"""20 pairs in the batch but only 19 unique ones: statistics are inferred, every orientation
	fails (fewer than 20 candidates), nothing is paired or rescued, and the output equals -P."""
	g, idx = ins_genome
	junk = C.PlantedGenome(998, [("x", 300)])
	pairs = long_inserts[:19] + [C.Pair("junk", C.Mate("junk", junk["x"][:150], "I" * 150, None, 0, 0, False, []),
		C.Mate("junk", junk["x"][150:], "I" * 150, None, 0, 0, False, []))]
	f1, f2 = _write(tmp_path, pairs)
	res = H.run_map(idx, [f1, f2])
	cands, pct, bnd = C.stderr_pestat(res.stderr)
	assert cands == [(0, 19, 0, 0)] and pct == [] and bnd == []
	assert not any(r.is_proper_pair for r in res.sam().records)
	unpaired = H.run_map(idx, [f1, f2], ["-P"])
	assert H.sam_body(res.stdout) == H.sam_body(unpaired.stdout)


@pytest.mark.parametrize("n_fr,rf_fails", [(400, False), (401, True)])
def test_min_orientation_ratio(ins_genome, tmp_path, n_fr, rf_fails):
	"""An orientation needs at least 5% of the largest orientation's candidates: 20 RF pairs
	next to 400 FR pairs are paired, next to 401 FR pairs they are not."""
	g, idx = ins_genome
	rng = H.Rng(79)
	pairs = []
	for i in range(n_fr):
		pairs.append(C.fr_pair(g, "f%d" % i, "chrA" if i < 300 else "chrB", 500 + (i % 300) * 450, rng.normal_int(300, 20), r1_rev=bool(i & 1)))
	rf = []
	for i in range(20):
		left = 1000 + i * 2900
		d = rng.normal_int(200, 15)
		rf.append(C.pair_at(g, "rf%d" % i, "chrB", left + 150, True, "chrB", left + 150 + d, False))
	f1, f2 = _write(tmp_path, pairs + rf)
	res = H.run_map(idx, [f1, f2])
	cands, pct, bnd = C.stderr_pestat(res.stderr)
	assert cands == [(0, n_fr, 20, 0)]
	assert [x[0] for x in bnd] == (["FR"] if rf_fails else ["FR", "RF"])
	by = _flags_by_pair(res.sam(), rf)
	assert all(a.is_proper_pair == b.is_proper_pair == (not rf_fails) for a, b in by.values())


def test_statistics_in_all_four_orientations(ins_genome, tmp_path):
	"""A library mixing all four orientations (30 pairs each from chrB: FR, FF with read 2 reversed, RR
	with read 1 reversed, RF with both) infers statistics for every orientation, each equal to the
	Python recomputation from the true 5' positions (helpers_c.reorient); FR and RF pairs within their
	bounds are proper. (Same-strand pairs are not checked for 0x2: see
	known-bugs/r2-same-strand-orientation-labels.md.)"""
	g, idx = ins_genome
	pairs = []
	for k, (orient, ins) in enumerate((("FR", 300), ("FF", 280), ("RR", 320), ("RF", 520))):
		for p in H.simulate_pairs(g, 30, seed=140 + k, insert_mean=ins, insert_sd=25, sub_rate=0.005, indel_rate=0.0,
				contigs=["chrB"], name_prefix=orient.lower()):
			pairs.append((p.name, orient) + C.reorient(p, orient))
	dists = [[], [], [], []]
	for name, orient, s1, s2, (code, dist) in pairs:
		dists[code].append(dist)
	assert [len(x) for x in dists] == [30, 30, 30, 30]
	ps = C.pestat(dists)
	f1 = H.write_fasta(str(tmp_path / "o1.fa"), [(x[0], x[2]) for x in pairs])
	f2 = H.write_fasta(str(tmp_path / "o2.fa"), [(x[0], x[3]) for x in pairs])
	res = H.run_map(idx, [f1, f2])
	cands, pct, bnd = C.stderr_pestat(res.stderr)
	assert cands == [(30, 30, 30, 30)]
	assert pct == [(C.DIR_NAMES[d], (s["p25"], s["p50"], s["p75"], "%.2f" % s["avg"], "%.2f" % s["std"])) for d, s in enumerate(ps)]
	assert bnd == [(C.DIR_NAMES[d], (s["lo"], s["hi"], "%.2f" % s["pct_lo"], "%.2f" % s["pct_hi"])) for d, s in enumerate(ps)]
	sam = res.sam()
	assert H.pair_invariant_errors(sam) == []
	prim = sam.primaries()
	for name, orient, s1, s2, (code, dist) in pairs:
		if orient in ("FR", "RF"):
			want = ps[code]["lo"] <= dist <= ps[code]["hi"]
			assert prim[(name, 1)].is_proper_pair == prim[(name, 2)].is_proper_pair == want, (name, dist)


# ---------------------------------------------------------------------------
# Statistics per batch
# ---------------------------------------------------------------------------


@pytest.mark.thread_invariance
def test_statistics_are_inferred_per_batch(ins_genome, tmp_path):
	"""Statistics are computed per batch. 70 pairs of ~300 bp then 20 of ~2000 bp: in one batch
	the 2000 bp pairs are outliers and not proper; with -K 21000 (70 pairs per batch) the second
	batch infers its own statistics and pairs them."""
	g, idx = ins_genome
	rng = H.Rng(80)
	a = [C.fr_pair(g, "a%d" % i, "chrA", 1000 + i * 1200, rng.normal_int(300, 10), r1_rev=bool(i & 1)) for i in range(70)]
	b = [C.fr_pair(g, "b%d" % i, "chrA", 90000 + i * 2400, rng.normal_int(2000, 10), r1_rev=bool(i & 1)) for i in range(20)]
	f1, f2 = _write(tmp_path, a + b)
	one = H.run_map(idx, [f1, f2])
	two = H.run_map(idx, [f1, f2], ["-K", "21000"])
	_, ps_all = _expected_stats(a + b)
	_, ps_a = _expected_stats(a)
	_, ps_b = _expected_stats(b)
	c1, _, bnd1 = C.stderr_pestat(one.stderr)
	c2, _, bnd2 = C.stderr_pestat(two.stderr)
	assert c1 == [(0, 90, 0, 0)] and c2 == [(0, 70, 0, 0), (0, 20, 0, 0)]  # one batch, then two
	assert [x[1][:2] for x in bnd1] == [(ps_all[1]["lo"], ps_all[1]["hi"])]
	assert [x[1][:2] for x in bnd2] == [(ps_a[1]["lo"], ps_a[1]["hi"]), (ps_b[1]["lo"], ps_b[1]["hi"])]
	by1, by2 = _flags_by_pair(one.sam(), a + b), _flags_by_pair(two.sam(), a + b)
	for p in a + b:
		want1 = ps_all[1]["lo"] <= p.dir_dist[1] <= ps_all[1]["hi"]
		assert by1[p.name][0].is_proper_pair == by1[p.name][1].is_proper_pair == want1, p
		assert by2[p.name][0].is_proper_pair and by2[p.name][1].is_proper_pair, p
	assert sum(by1[p.name][0].is_proper_pair for p in b) == 0


# ---------------------------------------------------------------------------
# -I
# ---------------------------------------------------------------------------

# option, (avg, std, hi, lo) by the documented defaults: std = round(avg/10), hi = avg + 4*std,
# lo = avg - 4*std raised to 1 (the code's max(); minibwa.1 says min(), see known-bugs)
INS_OPTS = [
	("500", (500, 50, 700, 300)),
	("500,20", (500, 20, 580, 420)),
	("500,50,620", (500, 50, 620, 300)),
	("500,50,650,380", (500, 50, 650, 380)),
	("333.4", (333, 33, 465, 201)),
	("100,50", (100, 50, 300, 1)),
]


@pytest.fixture(scope="module")
def ins_boundary_pairs(ins_genome):
	"""For every -I setting, FR pairs (both read orders) at lo-1, lo, hi and hi+1."""
	g, idx = ins_genome
	pairs, left = [], 1000
	seen = set()
	for opt, (avg, std, hi, lo) in INS_OPTS:
		for ins in (lo - 1, lo, hi, hi + 1):
			if ins < 1 or ins in seen:
				continue
			seen.add(ins)
			for r1_rev in (False, True):
				pairs.append(C.fr_pair(g, "i%d%s" % (ins, "r" if r1_rev else "f"), "chrA", left, ins, r1_rev=r1_rev))
				left += 1300
	return pairs


@pytest.mark.parametrize("opt,stats", INS_OPTS, ids=[o for o, _ in INS_OPTS])
def test_insert_size_option_bounds(ins_genome, ins_boundary_pairs, tmp_path, opt, stats):
	"""-I avg[,std[,max[,min]]] fixes the statistics (nothing is inferred although the batch has
	more than 20 pairs) and a pair is proper iff lo <= insert <= hi, with omitted values
	defaulting to std = round(avg/10), max = avg + 4*std and min = max(avg - 4*std, 1)."""
	g, idx = ins_genome
	avg, std, hi, lo = stats
	f1, f2 = _write(tmp_path, ins_boundary_pairs)
	res = H.run_map(idx, [f1, f2], ["-I", opt, "--rescue=0"])
	assert not C.has_pestat_report(res.stderr)
	sam = res.sam()
	assert H.pair_invariant_errors(sam) == []
	by = _flags_by_pair(sam, ins_boundary_pairs)
	n = 0
	for p in ins_boundary_pairs:
		a, b = by[p.name]
		dist = p.dir_dist[1]
		want = lo <= dist <= hi
		assert a.is_proper_pair == b.is_proper_pair == want, (opt, p, a.flag, b.flag)
		assert a.tlen == H.tlen_5prime(a, b)
		n += want
	assert n >= 4
