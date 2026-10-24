"""Area D: `minibwa mem`, the bwa-mem command line mapped onto minibwa's options.

From usage_mem() and main_mem(): mem always writes MD:Z; -W is the minimum
chaining score (map -m), -m the mate-rescue count (map --rescue), -S skips
rescue (--rescue 0), -T is map -s, -h is --xa, -z is --outs, -a outputs all
secondary records (--outn), -C copies comments (-y); -t -k -w -c -P -R -H -o
-5 -Y -I mean what they mean for map; options marked '*' are accepted and
ignored. So `mem ARGS` must write the same bytes as `map -b MD ARGS'`.
"""

import pytest

import helpers_d as D
import minibwa_helpers as H


class _Data(object):
	def __init__(self, **kw):
		self.__dict__.update(kw)


@pytest.fixture(scope="module")
def inputs(synth, tmp_path_factory):
	d = tmp_path_factory.mktemp("dmem")
	recs, info = D.mixed_reads(synth.genome, seed=47)
	se = D.write_records(str(d / "mixed.fq"), recs)
	pairs = H.simulate_pairs(synth.genome, 150, seed=121, insert_mean=300, insert_sd=30)
	p1, p2 = H.write_pairs(str(d / "p1.fq"), str(d / "p2.fq"), pairs)
	g = D.GradedRepeats()
	return _Data(se=se, p1=p1, p2=p2, graded=g, gidx=H.index_genome(g.genome, str(d)),
		gfq=H.write_fastq(str(d / "g.fq"), g.reads()), dir=d)


def _mem(index, reads, args, threads=None):
	return H.run_map(index, reads, list(args), threads=threads, subcommand="mem")


def test_mem_is_map_with_md(chrm):
	"""On the shipped chrM pairs, `mem` writes the same bytes as `map -b MD` (apart from @PG), and every
	mapped record carries MD:Z."""
	m = _mem(chrm.index, [chrm.read1, chrm.read2], [])
	assert H.strip_pg(m.stdout) == H.strip_pg(H.run_map(chrm.index, [chrm.read1, chrm.read2], ["-b", "MD"]).stdout)
	recs = m.sam().records
	assert all("MD" in r.tags for r in recs if not r.is_unmapped)
	assert m.sam().pg[0]["CL"].startswith("minibwa mem ")


@pytest.mark.parametrize("mem_args,map_args", [
	(["-t", "3", "-k", "21", "-w", "80", "-c", "100", "-W", "30", "-T", "40", "-m", "3", "-C", "-Y", "-5"],
		["-t", "3", "-k", "21", "-w", "80", "-c", "100", "-m", "30", "-s", "40", "--rescue", "3", "-y", "-Y", "-5"]),
	(["-S", "-P", "-I", "300,30"], ["--rescue", "0", "-P", "-I", "300,30"]),
	(["-m", "0", "-R", "@RG\\tID:m1", "-H", "@CO\\tmem"], ["--rescue", "0", "-R", "@RG\\tID:m1", "-H", "@CO\\tmem"]),
	(["-r", "1.5", "-y", "20", "-A", "3", "-B", "5", "-O", "7", "-E", "2", "-L", "9", "-U", "3", "-d", "50", "-x", "pacbio",
		"-p", "-j", "-q", "-K", "1000", "-V", "-M", "-u"], []),
], ids=["algorithm", "pairing", "rescue-rg-h", "ignored"])
def test_mem_options_map_onto_map(synth, inputs, mem_args, map_args):
	"""mem's options produce the output of the corresponding map options plus -b MD, on single-end reads
	(repeats, chimeras, comments) and on pairs; the options marked '*' change nothing."""
	for reads in (inputs.se, [inputs.p1, inputs.p2]):
		a = _mem(synth.index, reads, mem_args, threads=None if "-t" in mem_args else 2)
		b = H.run_map(synth.index, reads, map_args + ["-b", "MD"], threads=None if "-t" in map_args else 2)
		assert H.strip_pg(a.stdout) == H.strip_pg(b.stdout)


@pytest.mark.parametrize("mem_args,map_args", [(["-h", "3", "-z", "0.625"], ["--xa", "3", "--outs", "0.625"]),
	(["-h", "4", "-z", "0.625"], ["--xa", "4", "--outs", "0.625"]), (["-a", "-z", "0"], ["--outn", "1000000", "--outs", "0"]),
	(["-a", "-z", "0.75"], ["--outn", "1000000", "--outs", "0.75"])], ids=["h3", "h4", "a-z0", "a-z.75"])
def test_mem_xa_options(inputs, mem_args, map_args):
	"""-h sets the XA limit (--xa), -z the score ratio (--outs) and -a writes every secondary record that
	passes -z: on the graded repeats, `mem -a -z 0` writes all four secondaries of family A and all six of
	family B."""
	a = _mem(inputs.gidx, inputs.gfq, mem_args)
	b = H.run_map(inputs.gidx, inputs.gfq, map_args + ["-b", "MD"])
	assert H.strip_pg(a.stdout) == H.strip_pg(b.stdout)
	if mem_args == ["-a", "-z", "0"]:
		byq = a.sam().by_qname()
		assert [len(byq[q]) for q in ("q", "qr", "b")] == [5, 5, 7]


def test_mem_mask_level(tmp_path):
	"""-D sets the mask level: a 70 bp hit contained in the primary's query interval is a secondary (and
	then dropped as too weak) at the default -D 0.5 and at -D 0.9, but becomes a supplementary alignment
	when -D 1 lets a full overlap through."""
	rng = H.Rng(7)
	c = list(rng.dna(20000))
	read = list(c[5000:5150])
	read[70] = "ACGT"[("ACGT".index(read[70]) + 1) % 4]
	read = "".join(read)
	c[12000:12070] = list(read[:70])
	g = H.Genome([("c1", "".join(c))])
	idx = H.index_genome(g, str(tmp_path))
	fq = H.write_fastq(str(tmp_path / "d.fq"), [("d", read, "I" * 150)])
	for args in ([], ["-D", "0.9"]):
		recs = _mem(idx, fq, args).sam().records
		assert [(r.flag, r.pos, r.cigar) for r in recs] == [(0, 5001, "150M")]
	recs = _mem(idx, fq, ["-D", "1"]).sam().records
	assert [(r.flag, r.pos, r.cigar) for r in recs] == [(0, 5001, "150M"), (2048, 12001, "70M80H")]


def test_mem_verbosity(synth, inputs):
	"""-v 1 silences the [M::...] messages on stderr without changing the output; the default level writes
	them."""
	a = _mem(synth.index, inputs.se, [])
	b = _mem(synth.index, inputs.se, ["-v", "1"])
	assert any(l.startswith("[M::") for l in a.stderr.splitlines())
	assert not any(l.startswith("[M::") for l in b.stderr.splitlines())
	assert H.strip_pg(a.stdout) == H.strip_pg(b.stdout)


def test_mem_usage_and_errors(synth, inputs, tmp_path):
	"""`mem` without an index and reads prints its usage (with the adap defaults: -t 1, -k 19, -W 25,
	-m 10, -T 30, -h 5, -z 0.8, -D 0.5) to stderr and exits 1; an invalid -R exits 1; -o writes to a file."""
	u = H.run_minibwa(["mem"])
	assert u.returncode == 1 and u.stdout == ""
	assert u.stderr.startswith("Usage: minibwa mem [options] <in.idx> <in.fastq>\n")
	for s in ("    -t INT         number of threads [1]\n", "    -k INT         minimum seed length [19]\n",
			"    -W INT         min chaining score [25]\n", "    -m INT         mate rescue for up to INT candidates; 0 to skip rescue [10]\n",
			"    -T INT         suppress alignment with DP score lower than INT*{-A} [30]\n",
			"    -h INT         max secondary alignments at the XA tag [5]\n", "    -z FLOAT       skip XA if score < FLOAT*bestScore [0.8]\n",
			"    -D FLOAT       mask level [0.5]\n"):
		assert s in u.stderr, s
	r = H.run_minibwa(["mem", "-R", "@RG\\tSM:x", synth.index, inputs.se])
	assert r.returncode == 1 and r.stdout == "" and any(not l.startswith("[M::") for l in r.stderr.splitlines())
	out = str(tmp_path / "m.sam")
	r = _mem(synth.index, inputs.se, ["-o", out])
	with open(out) as f:
		assert r.stdout == "" and H.strip_pg(f.read()) == H.strip_pg(_mem(synth.index, inputs.se, []).stdout)
