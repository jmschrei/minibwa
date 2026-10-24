"""Area D: the PAF writer (mb_fmt_paf), -f/-a, and PAF with -u, -y, --outn/--outs and paired reads.

Each PAF line is derived independently from the SAM record of the same hit
(same options, same read set) and the genome: columns 1-9 and 12 from the
record's fields and clips, columns 10-11 (matching bases, alignment length
without ambiguous bases) recounted from the reference, and the tags.
"""

import struct

import pytest

import helpers_d as D
import minibwa_helpers as H


class _Data(object):
	def __init__(self, **kw):
		self.__dict__.update(kw)


@pytest.fixture(scope="module")
def mixed(synth, tmp_path_factory):
	d = tmp_path_factory.mktemp("dpaf")
	recs, info = D.mixed_reads(synth.genome, seed=43)
	fq = D.write_records(str(d / "mixed.fq"), recs)
	args = ["--outn", "3", "-b", "MD", "-y"]
	return _Data(recs=recs, info=info, fq=fq, reads=dict((r[0], r) for r in recs),
		sam=H.run_map(synth.index, fq, args), paf=H.run_map(synth.index, fq, ["-f"] + args),
		paf_u=H.run_map(synth.index, fq, ["-f", "-u"]), paf0=H.run_map(synth.index, fq, ["-f"]))


def _mlen_blen(ref, pos0, seq, cigar):
	"""(matching bases, alignment columns without ambiguous bases, ambiguous columns) of an alignment."""
	mlen = blen = ambi = 0
	r, q = pos0, 0
	for n, op in H.parse_cigar(cigar):
		if op in "M=X":
			for k in range(n):
				a, b = ref[r + k].upper(), seq[q + k].upper()
				if a not in "ACGT" or b not in "ACGT":
					ambi += 1
				else:
					blen += 1
					mlen += a == b
			r += n
			q += n
		elif op == "I":
			ambi += sum(c.upper() not in "ACGT" for c in seq[q:q + n])
			blen += sum(c.upper() in "ACGT" for c in seq[q:q + n])
			q += n
		elif op == "D":
			ambi += sum(c.upper() not in "ACGT" for c in ref[r:r + n])
			blen += sum(c.upper() in "ACGT" for c in ref[r:r + n])
			r += n
		elif op == "S":
			q += n
	return mlen, blen, ambi


def _expected_paf(rec, read, genome):
	"""The PAF line minibwa should write for a SAM record (except s1/cm/s2, which SAM does not carry)."""
	name, seq, qual, comment = read
	qlen = len(seq)
	if rec.is_unmapped:
		return "%s\t%d\t*\t*\t*\t*\t*\t*\t*\t0\t0\t0" % (name, qlen), None
	qs, qe = D.query_interval(rec, qlen)
	ref = genome[rec.rname]
	s = H.revcomp(seq) if rec.is_reverse else seq
	core = "".join("%d%s" % (n, op) for n, op in H.parse_cigar(rec.cigar) if op not in "SH")
	mlen, blen, ambi = _mlen_blen(ref, rec.pos0, s[rec.left_clip:], core)
	assert rec.tags["NM"] == blen - mlen + ambi
	cols = [name, qlen, qs, qe, rec.strand, rec.rname, len(ref), rec.pos0, rec.ref_end0, mlen, blen, rec.mapq]
	tags = ["tp:A:%s" % ("S" if rec.is_secondary else "P")]
	tail = ["NM:i:%d" % rec.tags["NM"], "AS:i:%d" % rec.tags["AS"], "ms:i:%d" % rec.tags["ms"], "md:i:%d" % rec.tags["md"],
		"cg:Z:" + core]
	if "MD" in rec.tags:
		tail.append("MD:Z:" + rec.tags["MD"])
	if comment:
		tail.append(comment)
	return "\t".join(str(c) for c in cols), (tags, tail)


def test_paf_lines_match_sam_records(synth, mixed):
	"""For every SAM record (same options) the PAF line has: name, read length, 0-based query start/end on
	the read as sequenced, strand, contig, contig length, 0-based target start/end, matching bases and
	alignment length (both recounted from the genome; ambiguous bases excluded; NM = blen-mlen+n_ambi),
	MAPQ; tp:A:P for primary and supplementary pieces and tp:A:S for secondaries; s1/cm on every line and
	s2 on tp:A:P lines only; NM/AS/ms/md equal to SAM; cg:Z = the CIGAR without clips; MD:Z with -b MD;
	the comment with -y. Unmapped reads give '*' columns and 0 0 0."""
	sam = mixed.sam.sam().records
	paf = mixed.paf.stdout.splitlines()
	assert len(sam) == len(paf)
	n_s = 0
	for rec, line in zip(sam, paf):
		head, rest = _expected_paf(rec, mixed.reads[rec.qname], synth.genome)
		f = line.split("\t")
		assert "\t".join(f[:12]) == head, (line[:120], rec.line[:120])
		if rest is None:
			assert len(f) == 12
			continue
		tags, tail = rest
		assert f[12] == tags[0]
		assert f[13].startswith("s1:i:") and int(f[13][5:]) > 0 and f[14].startswith("cm:i:") and int(f[14][5:]) > 0
		k = 15
		if f[12] == "tp:A:P":
			assert f[15].startswith("s2:i:") and int(f[15][5:]) >= 0
			k = 16
		else:
			n_s += 1
		assert f[k:] == tail, (line, rec.line)
	assert n_s >= 10


def test_paf_unmapped_and_u(mixed):
	"""-u drops the unmapped lines of PAF output and nothing else; without -y, no comment is written; PAF
	output has no header lines."""
	a = mixed.paf0.stdout.splitlines(True)
	assert not any(l.startswith("@") for l in a)
	kept = [l for l in a if l.split("\t")[5] != "*"]
	assert len(kept) == len(a) - 3
	assert mixed.paf_u.stdout.splitlines(True) == kept
	for l in a:
		assert not any(t.startswith(("BC:Z:", "XX:i:")) for t in l.split("\t"))


def test_paf_supplementary_pieces_are_primary(mixed):
	"""Each part of a chimeric read is a tp:A:P line with s2:i, and the parts cover disjoint query intervals."""
	byname = {}
	for l in mixed.paf0.stdout.splitlines():
		f = l.split("\t")
		byname.setdefault(f[0], []).append(f)
	for name, d in mixed.info.items():
		if d["cat"] != "chimeric":
			continue
		lines = byname[name]
		assert len(lines) == 2 and all(f[12] == "tp:A:P" and f[15].startswith("s2:i:") for f in lines)
		iv = sorted((int(f[2]), int(f[3])) for f in lines)
		assert iv[0][1] <= iv[1][0] + 5


def test_a_and_f_order(synth, mixed):
	"""-f selects PAF and the hidden -a selects SAM; the later option wins."""
	sam = H.run_map(synth.index, mixed.fq, ["-f", "-a"])
	paf = H.run_map(synth.index, mixed.fq, ["-a", "-f"])
	assert H.strip_pg(sam.stdout) == H.strip_pg(H.run_map(synth.index, mixed.fq).stdout)
	assert paf.stdout == mixed.paf0.stdout


def test_paired_names(synth, tmp_path):
	"""Paired reads get /1 and /2 appended in PAF (once, also when the input names end in /1 and /2);
	single-end reads do not; SAM drops an input /1 /2 suffix."""
	pairs = H.simulate_pairs(synth.genome, 25, seed=81, unique=True)
	for suffix in (False, True):
		f1, f2 = H.write_pairs(str(tmp_path / ("r1%d.fq" % suffix)), str(tmp_path / ("r2%d.fq" % suffix)), pairs, suffix=suffix)
		paf = H.run_map(synth.index, [f1, f2], ["-f"]).paf()
		names = [p.qname for p in paf]
		assert set(names) == set("%s/%d" % (p.name, i) for p in pairs for i in (1, 2))
		sam = H.run_map(synth.index, [f1, f2]).sam()
		assert set(r.qname for r in sam.records) == set(p.name for p in pairs)
	se = H.write_fastq(str(tmp_path / "se.fq"), [p.r1 for p in pairs])
	assert set(p.qname for p in H.run_map(synth.index, se, ["-f"]).paf()) == set(p.name for p in pairs)


def test_u_with_pairs(synth, tmp_path):
	"""In paired-end SAM and PAF, -u drops the record of an unmapped mate and keeps its mapped mate."""
	pairs = H.simulate_pairs(synth.genome, 30, seed=83, unique=True)
	rng = H.Rng(5)
	r2 = [p.r2 for p in pairs]
	for i in range(0, 30, 3):  # every third read 2 is random sequence
		r2[i] = H.SimRead(r2[i].name, rng.dna(len(r2[i].seq)), r2[i].qual, r2[i].segments)
	f1 = H.write_fastq(str(tmp_path / "1.fq"), [p.r1 for p in pairs])
	f2 = H.write_fastq(str(tmp_path / "2.fq"), r2)
	for fmt in ([], ["-f"]):
		a = H.run_map(synth.index, [f1, f2], fmt + ["-I", "300,30"]).stdout
		b = H.run_map(synth.index, [f1, f2], fmt + ["-I", "300,30", "-u"]).stdout
		body = [l for l in H.sam_body(a).splitlines(True)]
		if fmt:
			kept = [l for l in body if l.split("\t")[5] != "*"]
		else:
			kept = [l for l in body if not int(l.split("\t")[1]) & 4]
		assert len(kept) == len(body) - 10
		assert H.sam_body(b).splitlines(True) == kept


def test_chain_only_outs_uses_chain_score(index_of, tmp_path):
	"""Without base alignment (--chain-only, PAF), --outn writes secondary chains whose chaining score s1 is
	at least --outs times the primary's (the bound included); the s2:i of the primary is the best
	secondary chain score."""
	g = D.GradedRepeats()
	idx = index_of(g.genome)
	fq = H.write_fastq(str(tmp_path / "g.fq"), g.reads())
	every = H.run_map(idx, fq, ["--chain-only", "--outn", "100", "--outs", "0"]).paf()
	by = {}
	for p in every:
		by.setdefault(p.qname, []).append(p)
	assert all(len(v) >= 4 for v in by.values())
	for p in every:
		assert "cg" not in p.tags and "NM" not in p.tags
	for outs in ("0.9", "0.75", "0.5"):
		f = struct.unpack("f", struct.pack("f", float(outs)))[0]
		got = H.run_map(idx, fq, ["--chain-only", "--outn", "100", "--outs", outs]).paf()
		gb = {}
		for p in got:
			gb.setdefault(p.qname, []).append(p)
		for q, v in by.items():
			pri = v[0]
			assert pri.is_primary and all(not x.is_primary for x in v[1:])
			assert pri.tags["s2"] == max(x.tags["s1"] for x in v[1:])
			want = [pri.line] + [x.line for x in v[1:] if x.tags["s1"] >= f * pri.tags["s1"]]
			assert [x.line for x in gb[q]] == want


@pytest.mark.pinned
def test_pin_paf(mixed):
	"""Pins the current output: PAF of the mixed read set with --outn 3 -b MD -y."""
	H.check_pinned("d_mixed.paf.gz", mixed.paf.stdout)
