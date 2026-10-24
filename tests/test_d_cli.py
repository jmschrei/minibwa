"""Area D: the command line: `minibwa` dispatch, `minibwa map` option parsing, usage, errors, -o, -R, -H.

Expected messages and exit codes come from main.c and map-main.c, the usage
text and the README. Exit codes not stated anywhere (an unreadable read
file, an unwritable output file) are not asserted; see known-bugs/.
"""

import gzip
import os

import pytest

import helpers_d as D
import minibwa_helpers as H

SIGABRT = -6


@pytest.fixture(scope="module")
def reads(synth, tmp_path_factory):
	d = tmp_path_factory.mktemp("dcli")
	recs = H.simulate_reads(synth.genome, 30, seed=91, unique=True)
	fq = H.write_fastq(str(d / "r.fq"), recs)
	pairs = H.simulate_pairs(synth.genome, 30, seed=92, unique=True)
	il = H.write_interleaved(str(d / "il.fq"), pairs)
	return {"fq": fq, "il": il, "dir": d, "reads": recs}


def _run(args):
	return H.run_minibwa(args)


def test_no_arguments_prints_short_usage():
	"""`minibwa` alone prints the short command list to stdout and exits 0."""
	r = _run([])
	assert r.returncode == 0 and r.stderr == ""
	assert r.stdout.startswith("Usage: minibwa <command> <arguments>\nCommands:\n")
	for cmd in ("index", "map", "mem", "version"):
		assert "\n  %s " % cmd in r.stdout
	assert "Debugging" not in r.stdout


def test_help_prints_long_usage():
	"""`minibwa --help` prints the long command list (with the debugging commands) to stdout and exits 0."""
	r = _run(["--help"])
	assert r.returncode == 0
	assert r.stdout.startswith("Usage: minibwa <command> <arguments>\n")
	for cmd in ("index", "map", "mem", "version", "fastmap", "getref", "bench", "--help"):
		assert "\n    %s " % cmd in r.stdout


@pytest.mark.parametrize("argv", [["version"], ["map", "--version"]])
def test_version(argv):
	"""`minibwa version` and `minibwa map --version` print MB_VERSION from minibwa.h and exit 0."""
	r = _run(argv)
	assert r.returncode == 0
	assert r.stdout == D.mb_version() + "\n"


def test_unknown_command():
	"""An unknown command prints `ERROR: unknown command '<cmd>'` to stderr and exits 1."""
	r = _run(["frobnicate", "x"])
	assert r.returncode == 1 and r.stdout == ""
	assert r.stderr.strip() == "ERROR: unknown command 'frobnicate'"


def test_map_help_and_usage():
	"""`map --help` prints the usage to stdout and exits 0; `map` with fewer than two positional arguments
	prints it to stderr and exits 1. The defaults in it are those of the adap preset: -t 1, --outs 0.8,
	--xa 5 with '80%', -K 100m,1g."""
	h = _run(["map", "--help"])
	assert h.returncode == 0 and h.stdout.startswith("Usage: minibwa map [options] <in.idx> <in.fastq>\n")
	for line in ("    -t INT           number of worker threads [1]\n",
			"    --outs=FLOAT     output a secondary hit if score at least FLOAT*bestScore [0.8]\n",
			"    --xa=NUM         if <=NUM hits with score >80% of the best hit, output them to XA [5]\n",
			"    -K NUM1[,NUM2]   process NUM1-NUM2 bp of query sequences in a batch [100m,1g]\n",
			"    -o FILE          output file name [stdout]\n", "    -u               don't output unmapped reads\n",
			"    --outn=NUM       output up to {NUM,-N} secondary alignments [0]\n"):
		assert line in h.stdout, line
	for argv in (["map"], ["map", "only-an-index"]):
		u = _run(argv)
		assert u.returncode == 1 and u.stdout == ""
		assert u.stderr.startswith("Usage: minibwa map [options] <in.idx> <in.fastq>\n")
		assert u.stderr.count("\n") == h.stdout.count("\n")


def test_presets_are_applied_before_other_options():
	"""-x is applied in a first pass, so options given before it still take effect (the usage text shows
	the values in force when --help is reached): -t 7 and --outs 0.5 survive a later -x lr."""
	r = _run(["map", "-t", "7", "--outs", "0.5", "-x", "lr", "--help"])
	assert r.returncode == 0
	assert "number of worker threads [7]\n" in r.stdout
	assert "output a secondary hit if score at least FLOAT*bestScore [0.5]\n" in r.stdout
	assert ">50% of the best hit" in r.stdout
	assert "    -w NUM           bandwidth [500]\n" in r.stdout  # from -x lr


@pytest.mark.parametrize("argv,msg", [
	(["--no-such-option"], '[ERROR] unknown option in "--no-such-option"'),
	(["-Q"], '[ERROR] unknown option in "-Q"'),
	(["-x", "nonexistent"], "[ERROR] unknown preset 'nonexistent'"),
], ids=["long", "short", "preset"])
def test_map_option_errors(synth, reads, argv, msg):
	"""An unknown option or preset exits 1 with an error message and writes nothing to stdout."""
	r = _run(["map"] + argv + [synth.index, reads["fq"]])
	assert r.returncode == 1 and r.stdout == ""
	assert msg in r.stderr.splitlines()


def test_map_missing_option_argument(synth, reads):
	"""An option missing its argument exits 1 with `[ERROR] missing option argument`."""
	for argv in (["map", synth.index, reads["fq"], "-t"], ["map", synth.index, reads["fq"], "--outn"]):
		r = _run(argv)
		assert r.returncode == 1 and r.stdout == ""
		assert "[ERROR] missing option argument" in r.stderr.splitlines()


@pytest.mark.parametrize("cmd", ["map", "mem"])
def test_missing_index_aborts(reads, tmp_path, cmd):
	"""A missing index aborts (SIGABRT) with `[E::main_<cmd>] failed to load the index. ABORT!` before any
	output; so does --meth on an index built without --meth."""
	r = H.run_minibwa([cmd, str(tmp_path / "nope"), reads["fq"]])
	assert r.returncode == SIGABRT and r.stdout == ""
	assert "[E::main_%s] failed to load the index. ABORT!" % cmd in r.stderr


def test_meth_needs_meth_index(synth, reads):
	"""--meth loads <prefix>.meth.mbw; on an index built without --meth the load fails and map aborts."""
	r = H.run_minibwa(["map", "--meth", synth.index, reads["fq"]])
	assert r.returncode == SIGABRT and r.stdout == ""
	assert "[E::main_map] failed to load the index. ABORT!" in r.stderr


def test_unreadable_reads(synth, tmp_path):
	"""A read file that cannot be opened prints `ERROR: failed to open file '<name>'` and writes only the
	SAM header (the exit status is not asserted; see known-bugs/map-missing-reads-exit-0.md)."""
	missing = str(tmp_path / "missing.fq")
	r = _run(["map", synth.index, missing])
	assert any(l.startswith("ERROR: failed to open file '%s'" % missing) for l in r.stderr.splitlines())
	assert r.stdout and all(l.startswith("@") for l in r.stdout.splitlines())


def test_main_reports_version_and_command(synth, reads):
	"""After a successful command with arguments, main() writes `[M::main] Version: <MB_VERSION>` and
	`[M::main] CMD: <argv>` to stderr. (The timing lines are not asserted.)"""
	r = H.run_map(synth.index, reads["fq"])
	lines = r.stderr.splitlines()
	assert "[M::main] Version: %s" % D.mb_version() in lines
	assert "[M::main] CMD: %s" % " ".join(r.argv) in lines


def test_output_file(synth, reads, tmp_path):
	"""-o FILE writes the output to FILE and nothing to stdout; -o - writes to stdout. The bytes equal those
	written to stdout without -o (apart from the @PG command line)."""
	base = H.run_map(synth.index, reads["fq"])
	out = str(tmp_path / "o.sam")
	r = H.run_map(synth.index, reads["fq"], ["-o", out])
	assert r.stdout == ""
	with open(out) as f:
		text = f.read()
	assert H.strip_pg(text) == H.strip_pg(base.stdout)
	assert "\tCL:minibwa map -o %s " % out in text
	# run in a scratch directory, so a binary that treats "-" as a file name writes there, not into the tree
	r = H.run_minibwa(["map", "-o", "-", "-t", H.test_threads(), synth.index, reads["fq"]], cwd=str(tmp_path), check=True)
	assert H.strip_pg(r.stdout) == H.strip_pg(base.stdout)
	assert not os.path.exists(str(tmp_path / "-"))
	paf = str(tmp_path / "o.paf")
	r = H.run_map(synth.index, reads["fq"], ["-f", "-o", paf])
	with open(paf) as f:
		assert r.stdout == "" and f.read() == H.run_map(synth.index, reads["fq"], ["-f"]).stdout


@pytest.mark.parametrize("opt,name,flag_effect", [("--pe", "pe", True), ("--kalloc", "kalloc", False), ("--adap", "adap", False),
	("--long", "long", False)])
def test_yes_or_no(synth, reads, opt, name, flag_effect):
	"""Options taking yes/no accept yes, y, no and n; any other value prints
	`[WARNING] option '--<name>' only accepts 'yes' or 'no'.` and leaves the setting as it was. --pe=no
	(or n) stops interleaved reads with equal names from being paired; --pe=yes (or y) is the default."""
	base = H.strip_pg(H.run_map(synth.index, reads["il"]).stdout)
	bad = H.run_map(synth.index, reads["il"], ["%s=maybe" % opt])
	assert ("option '--%s' only accepts 'yes' or 'no'." % name) in bad.stderr
	assert H.strip_pg(bad.stdout) == base
	if not flag_effect:
		return
	for v in ("yes", "y"):
		r = H.run_map(synth.index, reads["il"], ["%s=%s" % (opt, v)])
		assert "WARNING" not in r.stderr and H.strip_pg(r.stdout) == base
	no = [H.strip_pg(H.run_map(synth.index, reads["il"], ["%s=%s" % (opt, v)]).stdout) for v in ("no", "n")]
	assert no[0] == no[1] != base
	assert all(r.flag & 0x1 for r in H.parse_sam(base).records)
	assert not any(r.flag & 0xC1 for r in H.parse_sam(no[0]).records)


# ---------------------------------------------------------------------------
# -R
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("rg,msg", [
	("RG\\tID:x", "[ERROR] the read group line is not started with @RG"),
	("@RG\tID:x", "[ERROR] the read group line contained literal <tab> characters -- replace with escaped tabs: \\t"),
	("@RG\\tSM:x", "[ERROR] no ID within the read group line"),
	("@RG\\tID:" + "a" * 256, "[ERROR] @RG:ID is longer than 255 characters"),
	("@RG\\tID:" + "b" * 300 + "\\tSM:x", "[ERROR] @RG:ID is longer than 255 characters"),
], ids=["no-@RG", "literal-tab", "no-ID", "ID-256", "ID-300"])
def test_rg_errors(synth, reads, rg, msg):
	"""-R must start with @RG, contain no literal tab, have an ID: field and an ID shorter than 256
	characters; otherwise `map` exits 1 with an error message and writes nothing to stdout."""
	r = H.run_map(synth.index, reads["fq"], ["-R", rg], check=False)
	assert r.returncode == 1 and r.stdout == ""
	assert msg in r.stderr.splitlines()


@pytest.mark.parametrize("rg,line,rid", [
	("@RG\\tID:" + "c" * 255, "@RG\tID:" + "c" * 255, "c" * 255),
	("@RG\\tSM:x\\tID:grp\\\\2\\tLB:y", "@RG\tSM:x\tID:grp\\2\tLB:y", "grp\\2"),
	("@RG\\tID:last", "@RG\tID:last", "last"),
], ids=["ID-255", "escapes", "ID-at-end"])
def test_rg_line_and_tag(synth, reads, rg, line, rid):
	"""-R with \\t and \\\\ escapes: the header gets the unescaped line after the @SQ lines and before @PG,
	and every record gets RG:Z with the ID (up to the next tab; 255 characters allowed)."""
	sam = H.run_map(synth.index, reads["fq"], ["-R", rg]).sam()
	assert [h for h in sam.header if h.startswith("@RG")] == [line]
	codes = [h[:3] for h in sam.header]
	assert codes.index("@RG") == codes.index("@PG") - 1 and codes.index("@RG") > max(i for i, c in enumerate(codes) if c == "@SQ")
	assert len(sam.records) == 30 and all(r.tags["RG"] == rid for r in sam.records)


# ---------------------------------------------------------------------------
# -H
# ---------------------------------------------------------------------------

def test_h_string(synth, reads):
	"""-H '@...' inserts the line (\\t, \\n escapes expanded) after @PG."""
	sam = H.run_map(synth.index, reads["fq"], ["-H", "@CO\\tline one\\n@CO\\tline two"]).sam()
	assert sam.header[-3:] == [sam.header[-3], "@CO\tline one", "@CO\tline two"]
	assert sam.header[-3].startswith("@PG\t")
	short = H.run_map(synth.index, reads["fq"], ["-H", "@CO"]).sam()  # too short to be a header line
	assert short.header[-1].startswith("@PG\t")


@pytest.mark.parametrize("gz", [False, True])
def test_h_file(synth, reads, tmp_path, gz):
	"""-H FILE (plain or gzip) inserts the lines of FILE that look like header lines (start with '@', a
	two-letter code and a tab) after @PG, in file order, and skips the others; a missing file adds
	nothing."""
	text = "@CO\tfirst\nnot a header\n@C\tshort code\n@CO\n@XYZ\tlong code\n\n@PG\tID:other\tPN:tool\n"
	path = str(tmp_path / ("h.txt.gz" if gz else "h.txt"))
	if gz:
		with gzip.open(path, "wt") as f:
			f.write(text)
	else:
		with open(path, "w") as f:
			f.write(text)
	sam = H.run_map(synth.index, reads["fq"], ["-H", path]).sam()
	i = [k for k, h in enumerate(sam.header) if h.startswith("@PG\t")][0]
	assert sam.header[i + 1:] == ["@CO\tfirst", "@PG\tID:other\tPN:tool"]
	none = H.run_map(synth.index, reads["fq"], ["-H", str(tmp_path / "missing.txt")]).sam()
	assert none.header[-1].startswith("@PG\t")


def test_h_and_r_ignored_in_paf(synth, reads):
	"""PAF output (-f) has no header, so -H and -R add nothing to it."""
	a = H.run_map(synth.index, reads["fq"], ["-f"]).stdout
	b = H.run_map(synth.index, reads["fq"], ["-f", "-H", "@CO\\tx", "-R", "@RG\\tID:z"]).stdout
	assert a == b and not a.startswith("@")
