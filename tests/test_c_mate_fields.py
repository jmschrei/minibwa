"""Area C: the mate half of the SAM writer (mb_fmt_sam, get_sam_pri): flags
0x1/0x2/0x8/0x20/0x40/0x80, RNEXT, PNEXT, TLEN, MC:Z and MQ:i, placement of an
unmapped mate; the /1 /2 suffix of paired names in PAF; and the options that
turn pairing off (-P, --pe, -x lr) or trace it (--dbg-qname).

Expected values come from the SAM specification and the truth of hand-built
pairs: FLAG bits from the placement of each mate; RNEXT/PNEXT from the mate's
primary record ('=' when on the same contig); an unmapped mate placed at its
mapped mate's position; TLEN the 5'-to-5' length (minibwa's documented
choice, H.tlen_5prime) when both are on one contig and 0 otherwise; MC:Z and
MQ:i the mate primary's CIGAR and MAPQ.
"""

import pytest

import minibwa_helpers as H
import helpers_c as C


class Scenario(object):
	def __init__(self, pair, map1, map2, proper):
		"""map1/map2: expected primary (contig, ts, reverse, cigar) or None for unmapped."""
		self.pair, self.map1, self.map2, self.proper = pair, map1, map2, proper

	@property
	def name(self):
		return self.pair.name


def _random_mate(name, seed, n=150):
	return C.Mate(name, C.PlantedGenome(seed, [("x", n)])["x"], C._qual(name + str(seed), n), None, 0, 0, False, [])


def _exp(m, cigar=None):
	return (m.contig, m.ts, m.rev, cigar or "%dM" % len(m))


def _build():
	pg = C.PlantedGenome(601, [("chrA", 80000), ("chrB", 40000)])
	rep = pg.random(150)
	pg.plant("chrA", 10250, rep)  # next to read 1 of "rep2"
	pg.plant("chrA", 30000, rep)  # the same contig, far away
	g = pg.genome()
	sc = []

	def fr(name, left, ins, r1_rev=False, proper=True, contig="chrA", len2=150):
		p = C.fr_pair(g, name, contig, left, ins, r1_rev=r1_rev, len2=len2)
		sc.append(Scenario(p, _exp(p.r1), _exp(p.r2), proper))
		return p

	fr("fr_f", 2000, 400)
	fr("fr_r", 3000, 400, r1_rev=True)
	fr("overlap", 4000, 200)  # the mates overlap by 100 bp
	fr("dovetail", 5000, 100, r1_rev=True)  # each mate extends past the other's 5' end
	fr("zero5", 6000, 1, proper=False)  # both 5' ends at one base: TLEN 0; 1 < lo (50)
	fr("far", 7000, 5000, proper=False)  # beyond hi (800)
	fr("short2", 13000, 300, len2=100)
	for k, (r1_rev, r2_rev) in enumerate(((False, False), (True, True))):  # same strand: FF/RR fail
		p = C.pair_at(g, "same%d" % k, "chrA", 14000 + k * 1000 + (150 if r1_rev else 0), r1_rev,
			"chrA", 14000 + k * 1000 + 400 + (150 if r2_rev else 0), r2_rev)
		sc.append(Scenario(p, _exp(p.r1), _exp(p.r2), False))
	p = C.pair_at(g, "inter", "chrA", 16000, False, "chrB", 2000, True)
	sc.append(Scenario(p, _exp(p.r1), _exp(p.r2), False))
	for k, r1_rev in enumerate((False, True)):
		p = C.fr_pair(g, "un2_%d" % k, "chrA", 17000 + k * 1000, 400, r1_rev=r1_rev)
		p.r2 = _random_mate(p.name, 610 + k)
		sc.append(Scenario(p, _exp(p.r1), None, False))
	p = C.fr_pair(g, "un1", "chrA", 19000, 400)
	p.r1 = _random_mate(p.name, 612)
	sc.append(Scenario(p, None, _exp(p.r2), False))
	p = C.Pair("un12", _random_mate("un12", 613), _random_mate("un12", 614))
	sc.append(Scenario(p, None, None, False))
	# read 2 in an exact repeat: one copy within the bounds of read 1, one 20 kbp away
	p = C.pair_at(g, "rep2", "chrA", 10000, False, "chrA", 10400, True)
	sc.append(Scenario(p, _exp(p.r1), _exp(p.r2), True))
	# chimeric read 1: 50 bp from chrB then 100 bp next to read 2 (read 2 is 100 bp)
	for k, (c5, p5) in enumerate((("chrB", 5000), ("chrA", 60000))):
		s = 21000 + k * 1000
		name = "chim%d" % k
		seq = g[c5][p5:p5 + 50] + g["chrA"][s:s + 100]
		r1 = C.Mate(name, seq, C._qual(name, 150), "chrA", s, s + 100, False, [])
		r2 = C.mate_at(g, name, "chrA", s + 300, True, length=100)
		p = C.Pair(name, r1, r2)
		p.supp = (c5, p5)
		sc.append(Scenario(p, ("chrA", s, False, "50S100M"), _exp(r2), True))
	single = C.mate_at(g, "single", "chrB", 9000, True)
	return g, sc, single


@pytest.fixture(scope="module")
def mf(index_of, tmp_path_factory):
	g, sc, single = _build()
	idx = index_of(g)
	d = tmp_path_factory.mktemp("mf")
	items = [s.pair for s in sc]
	items.insert(5, single)
	il = C.write_interleaved_c(str(d / "il.fq"), items)
	il_sfx = C.write_interleaved_c(str(d / "il_sfx.fq"), items, suffix=True)
	res = H.run_map(idx, il, ["--outn=1"])
	return g, sc, single, idx, il, il_sfx, res


def expect_flag(seg, me, mate, proper):
	"""FLAG of a primary record from the SAM specification; me/mate: (contig, ts, rev, cigar) or None."""
	f = 0x1 | (0x40 if seg == 1 else 0x80)
	if me is None:
		f |= 0x4
	elif me[2]:
		f |= 0x10
	if mate is None:
		f |= 0x8
	elif mate[2]:
		f |= 0x20
	if proper:
		f |= 0x2
	return f


def test_primary_records_match_spec(mf):
	"""Each mate's primary record: placement, FLAG, RNAME/POS of an unmapped read (its mate's),
	RNEXT/PNEXT, TLEN, MC:Z/MQ:i, and SEQ/QUAL of an unmapped read as given (not reversed)."""
	g, sc, single, idx, il, il_sfx, res = mf
	sam = res.sam()
	prim = sam.primaries()
	for s in sc:
		recs = (prim[(s.name, 1)], prim[(s.name, 2)])
		exps = (s.map1, s.map2)
		mates = (s.pair.r1, s.pair.r2)
		for i in (0, 1):
			r, m, me, mate = recs[i], recs[1 - i], exps[i], exps[1 - i]
			assert r.flag == expect_flag(i + 1, me, mate, s.proper), (s.name, i + 1, r.line[:100])
			if me is not None:
				assert (r.rname, r.pos0, r.is_reverse, r.cigar) == me, (s.name, r.line[:100])
			if me is None and mate is None:
				assert (r.rname, r.pos, r.rnext, r.pnext, r.tlen) == ("*", 0, "*", 0, 0)
			elif me is None:
				assert (r.rname, r.pos0, r.rnext, r.pnext) == (mate[0], mate[1], "=", mate[1] + 1)
			elif mate is None:
				assert (r.rnext, r.pnext) == ("=", r.pos)
			else:
				assert r.rnext == ("=" if me[0] == mate[0] else mate[0]) and r.pnext == mate[1] + 1
			if me is not None and mate is not None and me[0] == mate[0]:
				assert r.tlen == H.tlen_5prime(r, m)
			else:
				assert r.tlen == 0
			if me is not None and mate is not None:
				assert (r.tags["MC"], r.tags["MQ"]) == (m.cigar, m.mapq)
			else:
				assert "MC" not in r.tags and "MQ" not in r.tags
			if me is None:
				assert (r.seq, r.qual, r.mapq, r.cigar) == (mates[i].seq, mates[i].qual, 0, "*")
				assert not r.tags


def test_tlen_values(mf):
	"""TLEN by hand for the special geometries: +-400 for FR pairs, +-200 for overlapping mates,
	+-100 for a dovetailed pair (5'-to-5', not leftmost-to-rightmost: 200), 0 when the two 5'
	ends are one base, +-5000 for an improper pair on one contig."""
	g, sc, single, idx, il, il_sfx, res = mf
	prim = res.sam().primaries()
	want = {"fr_f": 400, "fr_r": -400, "overlap": 200, "dovetail": -100, "zero5": 0, "far": 5000, "short2": 300}
	for name, t in want.items():
		assert (prim[(name, 1)].tlen, prim[(name, 2)].tlen) == (t, -t), name
	d1, d2 = prim[("dovetail", 1)], prim[("dovetail", 2)]
	assert H.tlen_outer(d1, d2) == 200


def test_all_records_mate_fields(mf):
	"""Every record, including supplementary and secondary ones, carries its mate primary's
	position, strand (0x20) and CIGAR (MC:Z); TLEN is computed from the record itself. 0x2 is on
	every primary and supplementary record of a proper pair and never on a secondary one."""
	g, sc, single, idx, il, il_sfx, res = mf
	sam = res.sam()
	assert H.pair_invariant_errors(sam) == []
	prim = sam.primaries()
	by = dict((s.name, s) for s in sc)
	seen_supp = seen_sec = 0
	for r in sam.records:
		if r.qname == "single":
			continue
		s = by[r.qname]
		m = prim[(r.qname, 3 - r.segment)]
		if not r.is_unmapped and not m.is_unmapped and r.rname == m.rname:
			assert r.tlen == H.tlen_5prime(r, m), r.line[:100]
		else:
			assert r.tlen == 0
		if r.is_secondary:
			seen_sec += 1
			assert not r.is_proper_pair and r.seq == "*"
		else:
			assert r.is_proper_pair == s.proper, r.line[:100]
		if r.is_supplementary:
			seen_supp += 1
			assert (r.rname, r.pos0, r.cigar) == (s.pair.supp[0], s.pair.supp[1], "50M100H")
	assert seen_supp == 2 and seen_sec == 1
	sec = [r for r in sam.records if r.is_secondary][0]
	assert (sec.qname, sec.rname, sec.pos0, sec.flag) == ("rep2", "chrA", 30000, 0x100 | 0x1 | 0x10 | 0x80)
	assert sec.tlen == -(30149 - 10000 + 1)
	# MC:Z of read 2 is read 1's clipped CIGAR, which needs read 1's length (150), not read 2's (100)
	assert prim[("chim0", 2)].tags["MC"] == "50S100M"


def test_single_read_between_pairs(mf):
	"""A read without a mate in interleaved input has no pair flags, RNEXT '*', PNEXT 0, TLEN 0
	and no MC/MQ."""
	g, sc, single, idx, il, il_sfx, res = mf
	recs = [r for r in res.sam().records if r.qname == "single"]
	assert len(recs) == 1
	r = recs[0]
	assert (r.flag, r.rname, r.pos0, r.rnext, r.pnext, r.tlen) == (0x10, "chrB", 9000, "*", 0, 0)
	assert "MC" not in r.tags and "MQ" not in r.tags


def test_paf_mate_suffix(mf):
	"""PAF names paired reads name/1 and name/2 (input /1 /2 suffixes are not doubled) and leaves
	a single read's name alone; unmapped mates get a line of '*'."""
	g, sc, single, idx, il, il_sfx, res = mf
	for path in (il, il_sfx):
		paf = H.run_map(idx, path, ["-f"]).paf()
		names = set(p.qname for p in paf)
		want = set("%s/%d" % (s.name, k) for s in sc for k in (1, 2)) | set(["single"])
		assert names == want, path
		un = [p for p in paf if p.qname == "un12/2"]
		assert len(un) == 1 and un[0].is_unmapped and un[0].qlen == 150
		for p in paf:
			if p.qname.endswith("/1") or p.qname.endswith("/2"):
				s = [x for x in sc if x.name == p.qname[:-2]][0]
				m = s.pair.r1 if p.qname.endswith("/1") else s.pair.r2
				assert p.qlen == len(m)


def _se_view(recs):
	"""Alignment fields; an unmapped read's RNAME/POS are its mate's in paired output, so they are left out."""
	return [("unmapped",) if r.is_unmapped else (r.rname, r.pos, r.mapq, r.cigar, r.flag & 0x914, r.tags.get("NM"),
		r.tags.get("AS"), r.tags.get("ms")) for r in recs]


def test_no_pairing_keeps_single_end_alignments(mf):
	"""-P keeps the pair flags and mate fields but sets no 0x2 and changes no alignment: each
	read's records equal those of mapping the same file with --pe=no, where adjacent reads are
	not grouped into pairs at all (no 0x1, 0x40, 0x80, 0x8, 0x20; RNEXT '*')."""
	g, sc, single, idx, il, il_sfx, res = mf
	p = H.run_map(idx, il, ["-P", "--outn=1"]).sam()
	se = H.run_map(idx, il, ["--pe=no", "--outn=1"]).sam()
	assert H.pair_invariant_errors(p) == []
	assert not any(r.is_proper_pair for r in p.records)
	assert not any(r.flag & 0xE9 or r.rnext != "*" or r.tlen != 0 for r in se.records)
	assert any(r.is_paired for r in p.records)
	pv, sv = {}, {}
	for r in p.records:
		pv.setdefault((r.qname, r.segment or 1), []).append(r)
	order = {}
	for r in se.records:  # in --pe=no output, the first record set of a name is read 1
		if r.is_secondary or r.is_supplementary:
			sv[(r.qname, order[r.qname])].append(r)
			continue
		order[r.qname] = order.get(r.qname, 0) + 1
		sv[(r.qname, order[r.qname])] = [r]
	assert set(pv) == set(sv)
	for k in pv:
		assert _se_view(pv[k]) == _se_view(sv[k]), k


def test_long_read_preset_disables_pairing(mf):
	"""-x lr turns paired-end mode off (minibwa.1): interleaved reads are mapped as single reads;
	--pe=yes after the preset turns it back on and the unique FR pairs become proper."""
	g, sc, single, idx, il, il_sfx, res = mf
	lr = H.run_map(idx, il, ["-x", "lr"]).sam()
	assert not any(r.flag & 0xEB for r in lr.records)
	lrpe = H.run_map(idx, il, ["-x", "lr", "--pe=yes"]).sam()
	prim = lrpe.primaries()
	for name in ("fr_f", "fr_r", "overlap"):
		assert prim[(name, 1)].is_proper_pair and prim[(name, 2)].is_proper_pair


def test_dbg_qname(mf, threads):
	"""--dbg-qname (hidden) writes 'QN<TAB>name' lines naming the reads it maps and 'QP<TAB>name<TAB>thread'
	lines naming pairs it pairs, with a thread index below -t; the alignments do not change. How often
	a name is written is not asserted (an internal call count)."""
	g, sc, single, idx, il, il_sfx, res = mf
	dbg = H.run_map(idx, il, ["--outn=1", "--dbg-qname"])
	assert H.sam_body(dbg.stdout) == H.sam_body(res.stdout)
	qn, qp = set(), set()
	for l in dbg.stderr.splitlines():
		f = l.split("\t")
		if f[0] == "QN":
			qn.add(f[1])
		elif f[0] == "QP":
			qp.add(f[1])
			assert 0 <= int(f[2]) < threads
	assert qn == set([s.name for s in sc] + ["single"])
	assert qp and qp <= set(s.name for s in sc)
