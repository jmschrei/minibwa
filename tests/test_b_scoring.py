"""Area B: scoring options and MAPQ on hand-built reads, with expected values
worked out by hand from the documented scores (match -A, mismatch -B, N -1,
gap min(O1+l*E1, O2+l*E2), end bonus 10 for short reads), the -s/-m
thresholds (adjusted to the read length as minibwa.1 describes), and the
bwa-mem/minimap2 MAPQ formulas. The genome (helpers_b.BGenome) has a
two-copy unit differing at one base, a 400 bp unit in three exact copies and
a fourth copy with one substitution, and a stretch of unique sequence."""

import pytest

import minibwa_helpers as H
import helpers_b as B
from helpers_b import bgen  # noqa: F401 (fixture)


def _map(bgen, tmp_path, reads, args=(), name="b.fa"):
	return B.map_reads(bgen.index, tmp_path, reads, args, name=name)


def _mut(seq, i):
	return seq[:i] + B.flip_base(seq[i]) + seq[i + 1:]


def _clean_deletion(g, start, flank, dlen):
	"""The first position p >= start where deleting g[p+flank:p+flank+dlen] leaves no base that
	could slide across the deletion, so the alignment and any split are exactly at the cut."""
	p = start
	while g[p + flank - 1] == g[p + flank + dlen - 1] or g[p + flank] == g[p + flank + dlen]:
		p += 1
	return p


def test_mismatch_penalty_and_end_bonus(bgen, tmp_path):
	"""A mismatch 4 bases from the read end. Extending through it gains 3*2 - B and the end
	bonus 10; the alignment is clipped only when that is negative.
	-B 8: 100M, AS 99*2-8 = 190; -B 4: AS 194; -B 20: -14+10 < 0, so 96M4S with AS 192;
	-B 16: exactly 0, and a single-mismatch extension that ties goes to the end: 100M, AS 182."""
	u = bgen.unique_start
	s = _mut(bgen["b1"][u:u + 100], 96)
	want = {(): ("100M", 190, 1), ("-B", "4"): ("100M", 194, 1), ("-B", "20"): ("96M4S", 192, 0), ("-B", "16"): ("100M", 182, 1)}
	for args, (cigar, score, nm) in want.items():
		rec = B.one_record(_map(bgen, tmp_path, [("a", s)], args))
		assert (rec.pos0, rec.cigar, rec.tags["AS"], rec.tags["NM"]) == (u, cigar, score, nm), (args, rec.line[:200])
		assert rec.tags["ms"] == max(192, score)  # the best segment: the 96 matches or the whole read


def test_gap_open_and_extension(bgen, tmp_path):
	"""A 30 bp deletion between two 75 bp halves (150 matches, AS = 300 - gap cost):
	default dual affine min(12+2*30, 23+30) = 53 -> 247; -O 12 -E 2 (one value sets both
	pairs) 72 -> 228; -O 6,30 -E 3,1: min(96, 60) -> 240; -O 16 -E 2 -> 224.
	A 1 bp insertion (149 matches): 298 - 14 = 284; -O 6 -E 1: 291; -A 1: 149 - 14 = 135;
	-A 3: 447 - 14 = 433. ms charges any gap q+e in short-read mode: 300 - 14 = 286.
	(Pairs with O1 >= O2 or E1 <= E2 are not tested: base minibwa does not apply the
	documented cost to them, see known-bugs/r1-dual-gap-penalty-order.md.)"""
	g = bgen["b1"]
	u = _clean_deletion(g, bgen.unique_start + 1000, 75, 30)
	s = g[u:u + 75] + g[u + 105:u + 180]
	for args, score in (((), 247), (("-O", "12", "-E", "2"), 228), (("-O", "6,30", "-E", "3,1"), 240), (("-O", "16", "-E", "2"), 224)):
		rec = B.one_record(_map(bgen, tmp_path, [("d", s)], args))
		assert rec.pos0 == u and H.parse_cigar(rec.cigar) == [(75, "M"), (30, "D"), (75, "M")], (args, rec.line[:200])
		assert (rec.tags["AS"], rec.tags["NM"]) == (score, 30), args
	rec = B.one_record(_map(bgen, tmp_path, [("d", s)], ("-O", "12", "-E", "2")))
	assert rec.tags["ms"] == 300 - 14
	u = bgen.unique_start + 2000
	ins = [c for c in "ACGT" if c not in (g[u + 74], g[u + 75])][0]  # cannot slide left or right
	s = g[u:u + 75] + ins + g[u + 75:u + 149]
	for args, score in (((), 284), (("-O", "6", "-E", "1"), 291), (("-A", "1"), 135), (("-A", "3"), 433)):
		rec = B.one_record(_map(bgen, tmp_path, [("i", s)], args))
		assert (rec.pos0, rec.cigar, rec.tags["AS"], rec.tags["NM"]) == (u, "75M1I74M", score, 1), (args, rec.line[:200])


def test_match_and_mismatch_scores(bgen, tmp_path):
	"""Two mismatches in a 150 bp read: AS = 148*A - 2*B for -A 2 -B 8 (280), -A 3 -B 5
	(434) and -A 1 -B 1 (146)."""
	u = bgen.unique_start + 3000
	s = _mut(_mut(bgen["b1"][u:u + 150], 40), 110)
	for args, score in (((), 280), (("-A", "3", "-B", "5"), 434), (("-A", "1", "-B", "1"), 146)):
		rec = B.one_record(_map(bgen, tmp_path, [("m", s)], args))
		assert (rec.pos0, rec.cigar, rec.tags["AS"], rec.tags["NM"]) == (u, "150M", score, 2), (args, rec.line[:200])


def test_min_dp_score(bgen, tmp_path):
	"""50 matching bases then 100 unrelated ones: DP max 50*A. -s INT suppresses alignments
	below INT*A: kept at -s 50, unmapped at -s 51; with -A 3 the threshold scales (150 vs 153)."""
	g = bgen["b1"]
	u = bgen.unique_start + 4000
	junk = H.Rng(77).dna(100)
	junk = B.flip_base(g[u + 50]) + junk[1:]
	s = g[u:u + 50] + junk
	for args, mapped in ((("-s", "50"), True), (("-s", "51"), False), (("-A", "3", "-s", "50"), True), (("-A", "3", "-s", "51"), False)):
		rec = B.one_record(_map(bgen, tmp_path, [("s", s)], args))
		assert (not rec.is_unmapped) == mapped, (args, rec.line[:200])
		if mapped:
			assert rec.cigar == "50M100S" and rec.tags["AS"] == 50 * (3 if "-A" in args else 2)


def test_min_chain_score_and_mapq(bgen, tmp_path):
	"""A 30 bp exact match inside 150 bp: chain score 30, DP max 60. -m is adjusted for 150 bp
	to 40 - (40 - m)*2^(-50/1900): -m 31 -> 31 > 30, unmapped; -m 30 -> 30, mapped.
	MAPQ without a second hit: 40 * (1 - max(subsc, m')/30) * ln(60/2): -m 20 (m' 20) -> 45,
	default (25) -> 22, -m 30 -> 0, raised to 1 because DP max exceeds the (zero) second best."""
	g = bgen["b1"]
	u = bgen.unique_start + 5000
	rng = H.Rng(78)
	left, right = rng.dna(60), rng.dna(60)
	left = left[:-1] + B.flip_base(g[u - 1])
	right = B.flip_base(g[u + 30]) + right[1:]
	s = left + g[u:u + 30] + right
	for args, mapq in ((("-m", "20"), 45), ((), 22), (("-m", "30"), 1), (("-m", "31"), None)):
		rec = B.one_record(_map(bgen, tmp_path, [("c", s)], args))
		if mapq is None:
			assert rec.is_unmapped, (args, rec.line[:200])
			continue
		assert (rec.pos0, rec.cigar, rec.tags["AS"], rec.mapq) == (u, "60S30M60S", 60, mapq), (args, rec.line[:200])
	assert B.adap_value(31, 40, 150) == 31 and B.adap_value(30, 40, 150) == 30 and B.adap_value(25, 40, 150) == 25


def test_bandwidth(bgen, tmp_path):
	"""A 40 bp deletion in a 150 bp read. -w is adjusted for 150 bp (-w 20 -> 29, -w 5 -> 14)
	and the extension band is about 1.5 times that: -w 20 (and the default) bridge the
	deletion, -w 5 cannot and splits the read into a primary and a supplementary record."""
	g = bgen["b1"]
	u = _clean_deletion(g, bgen.unique_start + 6000, 75, 40)
	s = g[u:u + 75] + g[u + 115:u + 190]
	for args in ((), ("-w", "20")):
		recs = _map(bgen, tmp_path, [("w", s)], args).sam().records
		pri = [r for r in recs if r.is_primary]
		assert len(pri) == 1 and (pri[0].pos0, pri[0].cigar, pri[0].tags["AS"], pri[0].tags["NM"]) == (u, "75M40D75M", 300 - 63, 40), (args, [r.line[:120] for r in recs])
		assert not any(r.is_supplementary for r in recs)
	recs = _map(bgen, tmp_path, [("w", s)], ("-w", "5")).sam().records
	assert sorted((r.flag & 0x900, r.pos0, r.cigar.replace("H", "S")) for r in recs) == sorted(
		[(0, u, "75M75S"), (0x800, u + 115, "75S75M")]) or sorted((r.flag & 0x900, r.pos0, r.cigar.replace("H", "S")) for r in recs) == sorted(
		[(0x800, u, "75M75S"), (0, u + 115, "75S75M")]), [r.line[:120] for r in recs]
	assert all("SA" in r.tags for r in recs)


@pytest.mark.parametrize("tail,cigar,score", [(20, "120M20S", 240), (40, "120M30D40M", 80 + 240 - 53)])
def test_deletion_near_the_end(bgen, tmp_path, tail, cigar, score):
	"""120 bases, a 30 bp deletion, then a short tail. Aligning a 20 bp tail gains 40 - 53 =
	-13, and the end bonus (10) does not make up for it, so the tail is clipped; a 40 bp tail
	gains 80 - 53 = 27 and is aligned through the deletion."""
	g = bgen["b1"]
	u = _clean_deletion(g, bgen.unique_start + 13000, 120, 30)
	rec = B.one_record(_map(bgen, tmp_path, [("t", g[u:u + 120] + g[u + 150:u + 150 + tail])]))
	assert (rec.pos0, rec.cigar, rec.tags["AS"]) == (u, cigar, score), rec.line[:200]


def test_base_tag_choice(bgen, tmp_path):
	"""-b takes cs, ds, MD or md (one tag); any other value warns and writes cs."""
	u = bgen.unique_start + 14000
	s = [("b", _mut(bgen["b1"][u:u + 150], 70))]
	for args, tag, value in ((("-b", "md"), "MD", "70%s79" % bgen["b1"][u + 70]), (("-b", "cs"), "cs", ":70*%s%s:79" % (bgen["b1"][u + 70].lower(), B.flip_base(bgen["b1"][u + 70]).lower())),
			(("-b", "ds"), "ds", ":70*%s%s:79" % (bgen["b1"][u + 70].lower(), B.flip_base(bgen["b1"][u + 70]).lower()))):
		rec = B.one_record(_map(bgen, tmp_path, s, args))
		assert rec.tags[tag] == value and not set(rec.tags) & (set(["MD", "cs", "ds"]) - set([tag])), (args, rec.line[:200])
	res = _map(bgen, tmp_path, s, ("-b", "xyz"))
	assert any("-b" in l for l in B.warnings(res.stderr))  # a warning naming -b; wording not checked
	assert "cs" in B.one_record(res).tags


@pytest.mark.parametrize("length,mapq", [(60, 13), (100, 10), (150, 8), (250, 6)])
def test_mapq_short_two_copies(bgen, tmp_path, length, mapq):
	"""A read centred on the one base that tells two copies apart. Short-read (bwa-mem) MAPQ:
	x = 3/ln(blen), 6.02*x^2*(dp_max - dp_max2)/A rounded, dp_max - dp_max2 = 2+8 = 10, minus
	(int)(4.343*ln 2 + .499) = 3 for the near-equal hit: 60 bp 16-3, 100 bp 13-3, 150 bp 11-3,
	250 bp 9-3. md:i is 10 and XA lists the other copy with NM 1."""
	a, a2 = bgen.pos["A"], bgen.pos["A2"]
	lo = 500 - length // 2
	for name, unit, here, there in (("ra", bgen.unit_a, a, a2), ("rb", bgen.unit_a2, a2, a)):
		rec = B.one_record(_map(bgen, tmp_path, [(name, unit[lo:lo + length])]))
		assert (rec.pos0, rec.mapq, rec.tags["md"], rec.tags["n2"]) == (here + lo, mapq, 10, 1), rec.line[:200]
		assert rec.tags["XA"] == "b1,+%d,%dM,1;" % (there + lo + 1, length)


@pytest.mark.parametrize("length,args,mapq,md", [(400, (), 38, 52), (500, (), 33, 52), (600, ("-l", "1000"), 4, 10)])
def test_mapq_long_two_copies(bgen, tmp_path, length, args, mapq, md):
	"""Reads longer than -l are long reads: DP maxes are rescored as a*(mlen - b2*mismatches)
	with b2 = 25 (identity 1), so the other copy's 2*(L-1) - 2 becomes 2*(L-1-25) and md:i is 52;
	MAPQ blends the minimap2 and bwa-mem values: 400 bp: lr 30, sr 39 -> 30 + 9*2^(1-400/325)
	= 37.7 -> 38; 500 bp: lr 25, sr 36 -> 33. With -l 1000 a 600 bp read stays a short read:
	md 10 and bwa-mem MAPQ 7 - 3 = 4."""
	lo = 500 - length // 2
	rec = B.one_record(_map(bgen, tmp_path, [("lr", bgen.unit_a[lo:lo + length])], args))
	assert (rec.pos0, rec.cigar, rec.mapq, rec.tags["md"]) == (bgen.pos["A"] + lo, "%dM" % length, mapq, md), rec.line[:200]


def _r_read(bgen):
	return bgen.unit_r[125:275]  # the R4 copy differs at read offset 75


def _r_copies(bgen):
	"""(strand, 1-based position) of the read on each copy of R."""
	p = bgen.pos
	return {"R1": ("+", p["R1"] + 126), "R2": ("+", p["R2"] + 126), "R3": ("-", p["R3"] + 400 - 275 + 1), "R4": ("+", p["R4"] + 126)}


def test_secondaries_and_xa(bgen, tmp_path):
	"""A read from the 400 bp unit present three times exactly (one copy reverse-complemented)
	and once with one substitution: MAPQ 0, AS 300; n2:i:3 and XA lists the three other copies
	(NM 0, 0 and 1 for the R4 copy); with --outn the secondaries are SAM records (flag 256,
	SEQ '*', MAPQ 0) at those positions."""
	copies = _r_copies(bgen)
	res = _map(bgen, tmp_path, [("r", _r_read(bgen))], ("--outn", "10"))
	recs = res.sam().records
	pri = [r for r in recs if r.is_primary]
	assert len(pri) == 1 and pri[0].mapq == 0 and pri[0].tags["AS"] == 300 and pri[0].tags["n2"] == 3
	here = (pri[0].strand, pri[0].pos)
	assert here in [copies[k] for k in ("R1", "R2", "R3")]
	others = set(v for v in copies.values() if v != here)
	xa = set()
	for e in pri[0].tags["XA"].rstrip(";").split(";"):
		ctg, pos, cigar, nm = e.split(",")
		xa.add(((pos[0], int(pos[1:])), cigar, int(nm)))
	assert xa == set((v, "150M", 1 if v == copies["R4"] else 0) for v in others)
	sec = [r for r in recs if r.is_secondary]
	assert sorted((r.strand, r.pos) for r in sec) == sorted(others)
	assert all(r.mapq == 0 and r.seq == "*" and r.flag & 0x100 for r in sec)
	assert [r for r in sec if (r.strand, r.pos) == copies["R4"]][0].tags["AS"] == 290


def test_best_n_and_pri_ratio(bgen, tmp_path):
	"""-N caps the secondaries: -N 1 keeps one (n2 1), -N 0 none, and without a second hit
	the MAPQ formula gives 40*(1 - 150/150)*... = 0, raised to 1. -p 0.95 drops the R4 copy,
	whose chain is only the 75 bp seed beside the mismatch (75 < 0.95*150 and 75 + 2*19 < 150),
	and keeps the exact copies (n2 2)."""
	s = _r_read(bgen)
	rec = B.one_record(_map(bgen, tmp_path, [("r", s)], ("-N", "1")))
	assert rec.tags["n2"] == 1 and rec.tags["XA"].count(";") == 1 and rec.mapq == 0
	rec = B.one_record(_map(bgen, tmp_path, [("r", s)], ("-N", "0")))
	assert "n2" not in rec.tags and "XA" not in rec.tags and rec.mapq == 1 and rec.tags["md"] == 300
	rec = B.one_record(_map(bgen, tmp_path, [("r", s)], ("-p", "0.95")))
	assert rec.tags["n2"] == 2 and ",1;" not in rec.tags["XA"]
	rec = B.one_record(_map(bgen, tmp_path, [("r", s)], ("-p", "0")))
	assert rec.tags["n2"] == 3


def test_equal_hits_broken_by_read_name(bgen, tmp_path):
	"""Equal hits are ordered by a hash of the read name: the same sequence under 30 names
	lands on more than one of the exact copies, always with MAPQ 0 and AS 300, and each name's
	placement does not depend on the other reads or their order."""
	s = _r_read(bgen)
	copies = _r_copies(bgen)
	names = ["q%d" % i for i in range(30)]
	res1 = _map(bgen, tmp_path, [(n, s) for n in names], name="a.fa")
	res2 = _map(bgen, tmp_path, [(n, s) for n in reversed(names)] + [("extra", s)], name="b.fa")
	p1 = dict((r.qname, (r.strand, r.pos)) for r in res1.sam().records)
	p2 = dict((r.qname, (r.strand, r.pos)) for r in res2.sam().records)
	assert all(p1[n] == p2[n] for n in names)
	assert set(p1.values()) <= set(copies[k] for k in ("R1", "R2", "R3"))
	assert len(set(p1.values())) >= 2
	assert all(r.mapq == 0 and r.tags["AS"] == 300 for r in res1.sam().records)


def test_reference_n_run(bgen, tmp_path):
	"""A read across the 50 bp N run of b2 (read bases N where the reference has N): N counts
	-1 and as a mismatch in NM, as X with --eqx, as a match in cs and MD (both sides N),
	and not at all in PAF columns 10-11: AS 200-50 = 150, NM 50, mlen = blen = 100."""
	n = bgen.pos["N"]
	s = bgen["b2"][n - 50:n + 100]
	assert s[50:100] == "N" * 50
	rec = B.one_record(_map(bgen, tmp_path, [("n", s)], ("-b", "cs")))
	assert (rec.rname, rec.pos0, rec.cigar, rec.tags["AS"], rec.tags["NM"], rec.tags["cs"]) == ("b2", n - 50, "150M", 150, 50, ":150")
	rec = B.one_record(_map(bgen, tmp_path, [("n", s)], ("-b", "MD", "--eqx")))
	assert (rec.cigar, rec.tags["MD"]) == ("50=50X50=", "150")
	p = _map(bgen, tmp_path, [("n", s)], ("-f",)).paf()
	assert len(p) == 1 and (p[0].nmatch, p[0].alen, p[0].tags["NM"]) == (100, 100, 50)
