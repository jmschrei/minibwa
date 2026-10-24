"""Area C: the paired-end branch of the libminibwa API, mb_map_batch() with
MB_F_PE (reads 2i and 2i+1 are mates; the insert-size statistics come from
mb_opt_t), and the mb_hit_t fields proper_pair and rescued.

The C driver tests/c/drv_c_pe.c maps an interleaved file through the API and
writes each read's SAM-primary hit; the same reads mapped by `minibwa map -I`
with the same statistics must give the same placements, MAPQs and proper-pair
flags, since both run the same pairing on the same single-end hits.
"""

import os

import pytest

import minibwa_helpers as H
import helpers_c as C

PE = (400, 40, 560, 240)  # avg, std, hi, lo
UNSEED = C.unseedable_subs()


@pytest.fixture(scope="module")
def api_input(index_of, tmp_path_factory):
	pg = C.PlantedGenome(801, [("chrA", 60000), ("chrB", 20000)])
	g = pg.genome()
	pairs, kind = [], {}
	rng = H.Rng(81)
	for i in range(6):
		p = C.fr_pair(g, "u%d" % i, "chrA" if i < 4 else "chrB", 1000 + i * 3000, rng.normal_int(400, 30), r1_rev=bool(i & 1))
		pairs.append(p)
		kind[p.name] = "proper"
	for i in range(3):
		p = C.fr_pair(g, "s%d" % i, "chrA", 20000 + i * 3000, 400, r1_rev=bool(i & 1), subs2=UNSEED if i != 1 else (),
			subs1=UNSEED if i == 1 else ())
		pairs.append(p)
		kind[p.name] = "rescue%d" % (1 if i == 1 else 2)
	for name, ins in (("far", 700), ("far2", 580)):  # beyond hi (560); 580 is too close for rescue to clip
		p = C.fr_pair(g, name, "chrA", 32000 if name == "far" else 41000, ins)
		pairs.append(p)
		kind[p.name] = "improper"
	p = C.fr_pair(g, "rnd", "chrA", 35000, 400)
	p.r2 = C.Mate("rnd", C.PlantedGenome(802, [("x", 150)])["x"], "I" * 150, None, 0, 0, False, [])
	pairs.append(p)
	kind[p.name] = "unmapped2"
	odd = C.mate_at(g, "odd", "chrB", 15000, True)
	d = tmp_path_factory.mktemp("api")
	fa = str(d / "il.fa")
	recs = []
	for p in pairs:
		recs += [(p.name, p.r1.seq), (p.name, p.r2.seq)]
	recs.append(("odd", odd.seq))
	H.write_fasta(fa, recs)
	return g, index_of(g), pairs, kind, fa


def _driver(idx, fa, mode, out):
	H.assert_c_case("drv_c_pe", "pe_batch", idx, fa, PE[0], PE[1], PE[2], PE[3], mode, out)
	rows = []
	with open(out) as f:
		for line in f:
			x = line.rstrip("\n").split("\t")
			rows.append({"name": x[0], "i": int(x[1]), "n_hit": int(x[2]), "contig": x[3], "ts": int(x[4]),
				"te": int(x[5]), "strand": x[6], "mapq": int(x[7]), "proper": int(x[8]), "rescued": int(x[9]),
				"cigar": x[10], "n_res": int(x[13]), "n_pp": int(x[14])})
	return rows


def _cli_rows(res):
	"""SAM primaries in input order as the driver writes them."""
	out = []
	for r in res.sam().records:
		if not r.is_primary:
			continue
		core = "".join("%d%s" % (n, op) for n, op in r.cigar_ops if op not in "SH")
		if r.is_unmapped:
			out.append((r.qname, "*", -1, "*", 0, 0, "*"))
		else:
			out.append((r.qname, r.rname, r.pos0, r.strand, r.mapq, int(r.is_proper_pair), core))
	return out


def _api_view(rows):
	return [(x["name"], x["contig"], x["ts"], x["strand"], x["mapq"], x["proper"], x["cigar"]) for x in rows]


def test_api_pairing_matches_cli(api_input, tmp_path):
	"""mb_map_batch with MB_F_PE and pe_avg/std/hi/lo set gives every read the same primary
	placement, MAPQ and proper-pair status as `minibwa map -I avg,std,hi,lo`; the trailing odd
	read is mapped unpaired."""
	g, idx, pairs, kind, fa = api_input
	rows = _driver(idx, fa, "pe", str(tmp_path / "pe.tsv"))
	cli = H.run_map(idx, fa, ["-I", "%d,%d,%d,%d" % PE])
	assert _api_view(rows) == _cli_rows(cli)
	assert rows[-1]["name"] == "odd" and rows[-1]["proper"] == 0 and rows[-1]["n_hit"] == 1


def test_api_proper_pair_and_rescued_fields(api_input, tmp_path):
	"""proper_pair is set on both mates' primaries for pairs within [lo, hi] and on rescued mates;
	rescued is set exactly on the hits that mate rescue created (the unseedable mates) and the
	rescued mate lies at its true position."""
	g, idx, pairs, kind, fa = api_input
	rows = _driver(idx, fa, "pe", str(tmp_path / "pe.tsv"))
	by = {}
	for x in rows:
		by.setdefault(x["name"], []).append(x)
	truth = dict((p.name, p) for p in pairs)
	for name, rr in by.items():
		if name == "odd":
			continue
		k = kind[name]
		p = truth[name]
		assert [x["proper"] for x in rr] == ([1, 1] if k in ("proper", "rescue1", "rescue2") else [0, 0]), name
		for seg, (x, m) in enumerate(zip(rr, (p.r1, p.r2))):
			if k == "rescue%d" % (seg + 1):
				assert (x["rescued"], x["n_res"], x["n_hit"]) == (1, 1, 1), (name, x)
				assert (x["contig"], x["ts"], x["te"], x["strand"]) == (m.contig, m.ts, m.te, m.strand)
			else:
				assert x["n_res"] == 0, (name, x)
	assert sum(x["rescued"] for x in rows) == 3


@pytest.mark.parametrize("mode", ["nopair", "se"])
def test_api_without_pairing(api_input, tmp_path, mode):
	"""With MB_F_NO_PAIRING, or without MB_F_PE, no hit is proper or rescued and every read keeps
	its single-end placement and MAPQ, as `minibwa map -P` reports them; unseedable mates stay
	unmapped."""
	g, idx, pairs, kind, fa = api_input
	rows = _driver(idx, fa, mode, str(tmp_path / "np.tsv"))
	cli = H.run_map(idx, fa, ["-P"])
	assert _api_view(rows) == _cli_rows(cli)
	assert sum(x["n_hit"] == 0 for x in rows) == 4  # three unseedable mates and the random one
