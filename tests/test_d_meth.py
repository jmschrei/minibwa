"""Area D: --meth, mapping directional bisulfite reads against an `index --meth` index.

Per the README and minibwa.1: an unmethylated C on read 1 reads as T and an
unmethylated G on read 2 as A. Simulated directional pairs (original top
and original bottom strands, methylated CpGs kept as C) are mapped and
checked against their truth. In SAM orientation an OT read shows C->T
against the reference and an OB read G->A. Conversions score as matches
(AS) and are not counted in NM, while MD:Z and cs:Z show them as
mismatches; other differences count as usual.
"""

import os

import pytest

import helpers_d as D
import minibwa_helpers as H


class _Data(object):
	def __init__(self, **kw):
		self.__dict__.update(kw)


@pytest.fixture(scope="module")
def bs(synth, index_of, tmp_path_factory):
	d = tmp_path_factory.mktemp("dmeth")
	clean = D.bisulfite_pairs(synth.genome, 40, seed=111)
	noisy = [("e" + n, r1, r2, i) for n, r1, r2, i in D.bisulfite_pairs(synth.genome, 20, seed=112, sub_rate=0.02)]
	pairs = clean + noisy
	q = lambda s: "".join(chr(40 + (k * 3) % 30) for k in range(len(s)))
	r1 = H.write_fastq(str(d / "r1.fq"), [(n, a, q(a)) for n, a, b, i in pairs])
	r2 = H.write_fastq(str(d / "r2.fq"), [(n, b, q(b)) for n, a, b, i in pairs])
	idx = index_of(synth.genome, meth=True)
	return _Data(pairs=pairs, info=dict((p[0], p[3]) for p in pairs), r1=r1, r2=r2, index=idx,
		se=H.run_map(idx, r1, ["--meth", "-b", "MD"]), pe=H.run_map(idx, [r1, r2], ["--meth", "-b", "cs"]))


def _meth_of(rec, info, segment):
	"""Conversion visible in SAM orientation: OT pairs C->T, OB pairs G->A, whichever mate."""
	return "C2T" if info["strand"] == "OT" else "G2A"


def _truth_pos(info, segment, read_len=150):
	"""(strand, 0-based leftmost position) of a mate: the OT read 1 / OB read 2 is forward at the fragment start."""
	forward = (info["strand"] == "OT") == (segment in (0, 1))
	return ("+", info["frag_start"]) if forward else ("-", info["frag_end"] - read_len)


def _check_record(rec, genome, info, segment):
	strand, pos = _truth_pos(info, segment)
	assert (rec.rname, rec.strand) == (info["contig"], strand), rec.line[:120]
	assert abs(rec.unclipped_start0 - pos) <= 3, (rec.line[:120], info)
	meth = _meth_of(rec, info, segment)
	ref = genome[rec.rname]
	seq = rec.seq[rec.left_clip:len(rec.seq) - rec.right_clip]
	nm = H.compute_nm(ref, rec.pos0, rec.seq, rec.cigar, meth=meth)
	assert rec.tags["NM"] == nm, (rec.line[:120], nm)
	raw = H.compute_nm(ref, rec.pos0, rec.seq, rec.cigar)
	if rec.cigar == "150M":
		AS, ms = D.ungapped_scores(ref, rec.pos0, seq, meth=meth)
		assert (rec.tags["AS"], rec.tags["ms"]) == (AS, ms), rec.line[:120]
	return nm, raw


def test_single_end_bs_reads(synth, bs):
	"""Single-end reads are read 1 (C->T): OT reads map forward and OB reads reverse at their true positions
	with MAPQ >= 30; NM excludes conversions (0 for reads without other differences, whose AS is then 300)
	while MD:Z shows them as mismatches (the recomputed MD, without the final 0 base minibwa omits)."""
	sam = bs.se.sam()
	prim = sam.primaries()
	assert len(prim) == len(bs.pairs)
	n_conv = 0
	for (name, seg), rec in prim.items():
		info = bs.info[name]
		nm, raw = _check_record(rec, synth.genome, info, 0)
		assert rec.mapq >= 30
		md = H.compute_md(synth.genome[rec.rname], rec.pos0, rec.seq, rec.cigar)
		assert rec.tags["MD"] in (md, H.md_from_minibwa(md)), rec.line
		if not name.startswith("e"):
			assert nm == 0 and rec.tags["AS"] == 300 and rec.cigar == "150M"
		n_conv += raw - nm
	assert n_conv > 1000


def test_paired_bs_reads(synth, bs):
	"""Pairs: read 1 is C->T and read 2 G->A in read orientation, so both mates of an OT pair show C->T and
	both of an OB pair G->A in SAM orientation; both mates land at their true positions as a proper pair;
	NM excludes conversions; cs:Z shows them as *ct / *ga substitutions."""
	sam = bs.pe.sam()
	assert H.pair_invariant_errors(sam) == []
	prim = sam.primaries()
	assert len(prim) == 2 * len(bs.pairs)
	n_ct = n_ga = 0
	for (name, seg), rec in prim.items():
		info = bs.info[name]
		_check_record(rec, synth.genome, info, seg)
		assert rec.is_proper_pair
		assert rec.tags["cs"] == H.compute_cs(synth.genome[rec.rname], rec.pos0, rec.seq, rec.cigar), rec.line
		n_ct += rec.tags["cs"].count("*ct")
		n_ga += rec.tags["cs"].count("*ga")
		if info["strand"] == "OT":
			assert "*ga" not in rec.tags["cs"] or name.startswith("e")
		else:
			assert "*ct" not in rec.tags["cs"] or name.startswith("e")
	assert n_ct > 500 and n_ga > 500


def test_noisy_reads_count_other_mismatches(synth, bs):
	"""Differences other than the strand's conversion (sequencing errors after conversion) count in NM."""
	extra = 0
	for rec in bs.pe.sam().records:
		if rec.qname.startswith("e") and rec.is_primary:
			extra += rec.tags["NM"]
	assert extra >= 60


def _anchor_blocks(stderr):
	"""--dbg-qname/--dbg-anchor output at -t 1: [(qname, [(qs, strand, ctg, ts, len)])] in mapping order."""
	out = []
	for l in stderr.splitlines():
		f = l.split("\t")
		if f[0] == "QN":
			out.append((f[1], []))
		elif f[0] == "AC":
			assert out and out[-1][0] == f[1]
			out[-1][1].append((int(f[2]), f[3], f[4], int(f[5]), int(f[6])))
	return out


def test_anchors_split_at_non_conversion_mismatches(synth, bs):
	"""Seeds are found on converted sequences, so an anchor can span a base where the read has C and the
	reference T; minibwa splits anchors there. Every anchor (--dbg-anchor) therefore matches the reference
	base by base, allowing only the conversion of its strand (ref C/read T, or ref G/read A in reference
	orientation; read 1 forward and read 2 reverse are C->T), and is at least -k (19) bp long."""
	noisy = [p for p in bs.pairs if p[0].startswith("e")]
	il = H.write_fasta(os.path.join(os.path.dirname(bs.r1), "noisy_il.fa"), [x for n, a, b, i in noisy for x in ((n, a), (n, b))])
	r = H.run_map(bs.index, il, ["--meth", "--dbg-qname", "--dbg-anchor"], threads=1)
	blocks = _anchor_blocks(r.stderr)
	reads = {}
	for n, a, b, i in noisy:
		reads.setdefault(n, []).extend([a, b])
	seen = {}
	n_anchor = n_split = 0
	for name, anchors in blocks:
		k = seen.get(name, 0)
		seen[name] = k + 1
		read = reads[name][k]
		for qs, strand, ctg, ts, ln in anchors:
			assert ln >= 19
			q = read[qs:qs + ln]
			if strand == "-":
				q = H.revcomp(q)
			c2t = (k == 0) == (strand == "+")
			ref = synth.genome[ctg][ts:ts + ln].upper()
			for x, y in zip(ref, q):
				ok = x == y or (c2t and x == "C" and y == "T") or (not c2t and x == "G" and y == "A")
				assert ok, (name, k, qs, strand, ctg, ts, ln)
			n_anchor += 1
		n_split += len(anchors) > 1
	assert sorted(seen.values()) == [2] * len(noisy)
	assert n_anchor > 60 and n_split > 10


@pytest.mark.pinned
def test_pin_meth(bs):
	"""Pins the current output: --meth paired-end SAM with -b cs."""
	H.check_pinned("d_meth.sam.gz", H.strip_pg(bs.pe.stdout))
