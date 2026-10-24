"""Area D: output invariance of `minibwa map` across -t (kt_for, kt_pipeline), -K, --kalloc and --mmap.

minibwa's output is documented to be identical for every thread count; the
batch size -K must not change single-end output or paired-end output with
fixed insert-size statistics (-I), since only inferred statistics are
computed per batch. The kalloc allocator and the memory-mapped index
loaders are alternative implementations of the same thing. Every
comparison is of the full output bytes (apart from the @PG line).
"""

import pytest

import minibwa_helpers as H


class _Data(object):
	def __init__(self, **kw):
		self.__dict__.update(kw)


@pytest.fixture(scope="module")
def inputs(synth, tmp_path_factory):
	d = tmp_path_factory.mktemp("dinv")
	se = H.simulate_reads(synth.genome, 1000, length=(40, 250), seed=101, sub_rate=0.01, indel_rate=0.001, n_rate=0.002, qual=None)
	se += H.simulate_long_reads(synth.genome, 20, seed=102, qual=None)
	fq = H.write_fasta(str(d / "se.fa.gz"), se)
	pairs = H.simulate_pairs(synth.genome, 600, seed=103, insert_mean=350, insert_sd=40)
	p1, p2 = H.write_pairs(str(d / "p1.fq"), str(d / "p2.fq"), pairs)
	return _Data(fq=fq, p1=p1, p2=p2)


def _out(index, reads, args, threads):
	return H.strip_pg(H.run_map(index, reads, list(args), threads=threads).stdout)


@pytest.fixture(scope="module")
def se_base(synth, inputs):
	return _out(synth.index, inputs.fq, ["-b", "cs", "--outn", "2"], 1)


@pytest.mark.thread_invariance
@pytest.mark.parametrize("threads,args", [(2, ()), (3, ()), (4, ()), (8, ()), (3, ("-K", "5k")), (4, ("-K", "20k,40k")), (1, ("-K", "1",)),
	(4, ("--kalloc=yes",)), (2, ("--kalloc=yes", "-K", "3k")), (4, ("--kalloc=no",)), (2, ("--mmap",)), (3, ("--mmap=lite",))],
	ids=["t2", "t3", "t4", "t8", "t3-K5k", "t4-K20k,40k", "t1-K1", "t4-kalloc", "t2-kalloc-K3k", "t4-nokalloc", "t2-mmap", "t3-mmap-lite"])
def test_single_end_invariance(synth, inputs, se_base, threads, args):
	"""Single-end SAM (short and long reads, --outn 2, -b cs) is byte-identical to -t 1 for -t 2/3/4/8, for
	batches of 1 bp to 40 kbp (-K), with kalloc on or off, and with the index memory-mapped (--mmap,
	--mmap=lite)."""
	assert _out(synth.index, inputs.fq, ["-b", "cs", "--outn", "2"] + list(args), threads) == se_base


@pytest.mark.thread_invariance
def test_paf_invariance(synth, inputs):
	"""PAF output does not depend on -t or -K either."""
	a = _out(synth.index, inputs.fq, ["-f"], 1)
	assert _out(synth.index, inputs.fq, ["-f"], 5) == a
	assert _out(synth.index, inputs.fq, ["-f", "-K", "7k"], 3) == a


@pytest.mark.thread_invariance
def test_paired_end_invariance(synth, inputs):
	"""Paired-end SAM with inferred insert sizes is identical for -t 1 and 4 at the default -K; with -I it is
	also identical across -K values and kalloc."""
	a = _out(synth.index, [inputs.p1, inputs.p2], [], 1)
	assert _out(synth.index, [inputs.p1, inputs.p2], [], 4) == a
	b = _out(synth.index, [inputs.p1, inputs.p2], ["-I", "350,40"], 1)
	for t, extra in ((4, ["-K", "6k"]), (2, ["-K", "20k", "--kalloc=yes"]), (3, ["--mmap"])):
		assert _out(synth.index, [inputs.p1, inputs.p2], ["-I", "350,40"] + extra, t) == b


@pytest.mark.thread_invariance
def test_chrm_invariance(chrm):
	"""The shipped chrM pairs (real reads, inferred insert sizes): -t 1 and -t 4 agree, in SAM and PAF."""
	for args in ([], ["-f"]):
		assert _out(chrm.index, [chrm.read1, chrm.read2], args, 4) == _out(chrm.index, [chrm.read1, chrm.read2], args, 1)


@pytest.mark.pinned
def test_pin_chrm(chrm):
	"""Pins the current output: SAM of the shipped chrM pairs with -b MD (the README's example)."""
	H.check_pinned("d_chrm.sam.gz", _out(chrm.index, [chrm.read1, chrm.read2], ["-b", "MD"], None))


@pytest.mark.pinned
def test_pin_single_end(se_base):
	"""Pins the current output: single-end SAM of 1,020 simulated short and long reads (FASTA) with -b cs --outn 2."""
	H.check_pinned("d_se.sam.gz", se_base)
