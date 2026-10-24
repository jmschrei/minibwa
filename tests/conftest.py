"""Session fixtures shared by every test file.

minibwa_bin     path of the binary under test (built out of tree unless $MINIBWA_BIN)
threads         $MINIBWA_TEST_THREADS (default 1); run_map() uses it unless told otherwise
index_dir       a session-lived directory for indexes (pass it to index_genome)
index_of        index_of(genome_or_fasta, meth=False) -> index prefix, cached
chrm            the shipped test data: chrM genome, its index and the two read files
synth           the standard synthetic genome (see make_synth) and its index
tree_before     build products present in the source tree before the session started
"""

import os

import pytest

import minibwa_build
import minibwa_helpers as H

_TREE_BEFORE = {}


def pytest_configure(config):
	global _TREE_BEFORE
	_TREE_BEFORE = minibwa_build.tree_products()


def pytest_report_header(config):
	lines = ["minibwa source: %s" % minibwa_build.src_dir()]
	try:
		lines.append("minibwa binary: %s" % minibwa_build.minibwa_bin())
	except minibwa_build.BuildError as e:
		lines.append("minibwa binary: BUILD FAILED (%s)" % str(e).splitlines()[0])
	lines.append("CC=%s  MINIBWA_TEST_THREADS=%d" % (minibwa_build.compiler(), H.test_threads()))
	return lines


@pytest.fixture(scope="session")
def minibwa_bin():
	try:
		return minibwa_build.minibwa_bin()
	except minibwa_build.BuildError as e:
		pytest.fail("cannot build minibwa: %s" % e, pytrace=False)


@pytest.fixture(scope="session")
def threads():
	return H.test_threads()


@pytest.fixture(scope="session")
def tree_before():
	return dict(_TREE_BEFORE)


@pytest.fixture(scope="session")
def index_dir(tmp_path_factory):
	return str(tmp_path_factory.mktemp("indexes"))


@pytest.fixture(scope="session")
def index_of(minibwa_bin, index_dir):
	def f(genome_or_fasta, meth=False, extra_args=()):
		return H.index_genome(genome_or_fasta, index_dir, meth=meth, binary=minibwa_bin, extra_args=extra_args)
	return f


class _Data(object):
	def __init__(self, **kw):
		self.__dict__.update(kw)


@pytest.fixture(scope="session")
def chrm(index_of):
	"""test/chrM-human.fa.gz (one 16,569 bp contig) and 1,000 read pairs of 151 bp."""
	fasta = os.path.join(H.DATA_DIR, "chrM-human.fa.gz")
	return _Data(fasta=fasta, genome=H.read_fasta(fasta), index=index_of(fasta),
		read1=os.path.join(H.DATA_DIR, "chrM-read_1.fa.gz"), read2=os.path.join(H.DATA_DIR, "chrM-read_2.fa.gz"))


def make_synth():
	"""The standard synthetic genome: two contigs (60 and 40 kbp) with two exact-repeat
	families of three 300-1000 bp copies, three low-complexity tracts, two N runs and
	two short contigs (20-200 bp). Changing it changes every pin built on it."""
	return H.make_genome(seed=11, contig_lengths=(60000, 40000), repeat_families=2, repeat_copies=3,
		low_complexity=3, n_runs=2, short_contigs=2)


@pytest.fixture(scope="session")
def synth(index_of, index_dir):
	g = make_synth()
	fasta = g.write_fasta(os.path.join(index_dir, "synth.fa"))
	return _Data(genome=g, fasta=fasta, index=index_of(g))
