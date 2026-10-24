"""Locate or build the minibwa binary, libminibwa.a and the C test drivers.

Nothing here writes to the source tree. Sources are copied to a cache
directory keyed by a hash of the sources and the compiler, and built there
with the Makefile's own flags. A file lock serializes builds so parallel
pytest workers and concurrent test runs do not race.

Environment variables:
	MINIBWA_BIN            use this binary instead of building one
	MINIBWA_TEST_SRC_DIR   source tree to build from (default: repository root)
	MINIBWA_TEST_CACHE_DIR build cache (default: $XDG_CACHE_HOME/minibwa-tests
	                       or ~/.cache/minibwa-tests)
	MINIBWA_TEST_MAKE_JOBS parallel jobs for make (default: 4)
	CC                     C compiler (default: gcc)
"""

import fcntl
import hashlib
import os
import shutil
import subprocess
import sys
import tempfile

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(TESTS_DIR)
C_DIR = os.path.join(TESTS_DIR, "c")

# Build products that `make` would leave in a source tree; used to check that
# the harness never builds in place.
TREE_PRODUCT_PATTERNS = ("*.o", "*.a", "minibwa", "a.out", "api-test/mbmap-one", "api-test/mbmap-batch")


class BuildError(RuntimeError):
	pass


def src_dir():
	return os.path.abspath(os.environ.get("MINIBWA_TEST_SRC_DIR") or REPO_ROOT)


def compiler():
	return os.environ.get("CC") or "gcc"


def make_jobs():
	try:
		return max(1, int(os.environ.get("MINIBWA_TEST_MAKE_JOBS", "4")))
	except ValueError:
		return 4


def cache_root():
	d = os.environ.get("MINIBWA_TEST_CACHE_DIR")
	if not d:
		base = os.environ.get("XDG_CACHE_HOME") or os.path.join(os.path.expanduser("~"), ".cache")
		d = os.path.join(base, "minibwa-tests")
	os.makedirs(d, exist_ok=True)
	return d


def source_files(root):
	"""Relative paths of the files `make` reads: Makefile, top-level *.c/*.h, mimalloc/."""
	out = []
	for name in sorted(os.listdir(root)):
		p = os.path.join(root, name)
		if os.path.isfile(p) and (name == "Makefile" or name.endswith(".c") or name.endswith(".h")):
			out.append(name)
	mi = os.path.join(root, "mimalloc")
	for dirpath, dirnames, filenames in os.walk(mi):
		dirnames.sort()
		for name in sorted(filenames):
			out.append(os.path.relpath(os.path.join(dirpath, name), root))
	if "Makefile" not in out:
		raise BuildError("no Makefile in %s; set MINIBWA_TEST_SRC_DIR to a minibwa source tree" % root)
	return out


def _hash_files(h, root, rels):
	for rel in rels:
		h.update(rel.replace(os.sep, "/").encode() + b"\0")
		with open(os.path.join(root, rel), "rb") as f:
			h.update(hashlib.sha256(f.read()).digest())


def source_hash(root=None):
	root = root or src_dir()
	h = hashlib.sha256()
	_hash_files(h, root, source_files(root))
	return h.hexdigest()


def compiler_id(cc=None):
	"""A string that changes when the compiler does: its resolved path and version banner."""
	cc = cc or compiler()
	exe = cc.split()[0]
	path = shutil.which(exe) or exe
	try:
		ver = subprocess.run(cc.split() + ["--version"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
			timeout=60).stdout.decode("utf-8", "replace")
	except (OSError, subprocess.SubprocessError) as e:
		raise BuildError("cannot run compiler %r: %s" % (cc, e))
	return "%s\n%s\n%s" % (cc, os.path.realpath(path), ver)


def build_key(root=None, cc=None):
	h = hashlib.sha256()
	h.update(source_hash(root).encode())
	h.update(compiler_id(cc).encode())
	h.update(("%s-%s" % (sys.platform, os.uname().machine)).encode())
	return h.hexdigest()[:20]


class _Lock(object):
	def __init__(self, path):
		self.path = path
		self.fd = None

	def __enter__(self):
		self.fd = os.open(self.path, os.O_CREAT | os.O_RDWR, 0o644)
		fcntl.flock(self.fd, fcntl.LOCK_EX)
		return self

	def __exit__(self, *exc):
		fcntl.flock(self.fd, fcntl.LOCK_UN)
		os.close(self.fd)
		self.fd = None


def _run(cmd, cwd, log):
	p = subprocess.run(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=900)
	with open(log, "ab") as f:
		f.write(("$ %s\n" % " ".join(cmd)).encode())
		f.write(p.stdout)
	if p.returncode != 0:
		tail = p.stdout.decode("utf-8", "replace")[-4000:]
		raise BuildError("command failed (exit %d): %s\n%s" % (p.returncode, " ".join(cmd), tail))


def build_dir(root=None, cc=None):
	"""Build minibwa and libminibwa.a out of tree; return the build directory.

	The directory holds a copy of the sources plus `minibwa`, `libminibwa.a`
	and `build.log`. It is reused while the sources and the compiler stay the
	same.
	"""
	root = root or src_dir()
	cc = cc or compiler()
	key = build_key(root, cc)
	cache = cache_root()
	final = os.path.join(cache, "build-" + key)
	marker = os.path.join(final, ".complete")
	if os.path.exists(marker):
		return final
	with _Lock(os.path.join(cache, "build-" + key + ".lock")):
		if os.path.exists(marker):
			return final
		if os.path.exists(final):
			shutil.rmtree(final)
		tmp = tempfile.mkdtemp(prefix="tmp-build-" + key + "-", dir=cache)
		try:
			for rel in source_files(root):
				dst = os.path.join(tmp, rel)
				os.makedirs(os.path.dirname(dst), exist_ok=True)
				shutil.copy2(os.path.join(root, rel), dst)
			if build_key(tmp, cc) != key:
				raise BuildError("sources changed while being copied from %s; retry" % root)
			log = os.path.join(tmp, "build.log")
			_run(["make", "-j%d" % make_jobs(), "CC=" + cc, "minibwa", "libminibwa.a"], tmp, log)
			for prod in ("minibwa", "libminibwa.a"):
				if not os.path.exists(os.path.join(tmp, prod)):
					raise BuildError("make did not produce %s in %s" % (prod, tmp))
			with open(os.path.join(tmp, "build.info"), "w") as f:
				f.write("source: %s\nkey: %s\ncompiler: %s\n" % (root, key, compiler_id(cc)))
			open(os.path.join(tmp, ".complete"), "w").close()
			os.replace(tmp, final)
		except BaseException:
			shutil.rmtree(tmp, ignore_errors=True)
			raise
	return final


def minibwa_bin():
	"""Path of the binary under test: $MINIBWA_BIN, else an out-of-tree build."""
	b = os.environ.get("MINIBWA_BIN")
	if b:
		b = os.path.abspath(b)
		if not (os.path.isfile(b) and os.access(b, os.X_OK)):
			raise BuildError("MINIBWA_BIN=%s is not an executable file" % b)
		return b
	return os.path.join(build_dir(), "minibwa")


# Objects that `make minibwa` builds but does not put in libminibwa.a, which C drivers may
# call: bseq.o (mb_bseq_*) and kthread.o (kt_for, kt_pipeline). They are archived into
# libmbt_extra.a in the build directory, so a driver links only the members it needs and a
# driver that includes the .c file itself still links without duplicate symbols.
EXTRA_OBJS = ("bseq.o", "kthread.o")
EXTRA_LIB = "libmbt_extra.a"


def extra_lib(d=None):
	"""Path of libmbt_extra.a (EXTRA_OBJS from the out-of-tree build), created on first use."""
	d = d or build_dir()
	lib = os.path.join(d, EXTRA_LIB)
	if os.path.exists(lib):
		return lib
	with _Lock(lib + ".lock"):
		if os.path.exists(lib):
			return lib
		objs = [os.path.join(d, o) for o in EXTRA_OBJS]
		missing = [o for o in objs if not os.path.exists(o)]
		if missing:
			raise BuildError("the build in %s has no %s" % (d, ", ".join(missing)))
		tmp = lib + ".tmp%d" % os.getpid()
		cmd = (os.environ.get("AR") or "ar").split() + ["rcs", tmp] + objs
		p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=120)
		if p.returncode != 0:
			raise BuildError("command failed (exit %d): %s\n%s" % (p.returncode, " ".join(cmd), p.stdout.decode("utf-8", "replace")))
		os.replace(tmp, lib)
	return lib


def libminibwa():
	"""(path to libminibwa.a, include directory), both from the out-of-tree build."""
	d = build_dir()
	return os.path.join(d, "libminibwa.a"), d


def c_driver_sources(name):
	src = os.path.join(C_DIR, name + ".c")
	if not os.path.exists(src):
		raise BuildError("no C driver %s" % src)
	return [src, os.path.join(C_DIR, "minibwa_test.h")]


def c_driver(name, extra_cflags=()):
	"""Compile tests/c/<name>.c against libminibwa.a and libmbt_extra.a; return the executable's path.

	The driver sees minibwa's headers (minibwa.h and the internal ones) and
	tests/c/minibwa_test.h. It is rebuilt when the driver, the test header,
	the library sources or the compiler change.
	"""
	lib, inc = libminibwa()
	xlib = extra_lib(inc)
	cc = compiler()
	h = hashlib.sha256()
	h.update(os.path.basename(inc).encode())
	h.update(" ".join(extra_cflags).encode())
	h.update(EXTRA_LIB.encode())
	for p in c_driver_sources(name):
		with open(p, "rb") as f:
			h.update(hashlib.sha256(f.read()).digest())
	out_dir = os.path.join(inc, "drivers")
	os.makedirs(out_dir, exist_ok=True)
	exe = os.path.join(out_dir, "%s-%s" % (name, h.hexdigest()[:16]))
	if os.path.exists(exe):
		return exe
	with _Lock(exe + ".lock"):
		if os.path.exists(exe):
			return exe
		tmp = exe + ".tmp%d" % os.getpid()
		cmd = cc.split() + ["-std=gnu99", "-g", "-O2", "-Wall"] + list(extra_cflags) + \
			["-I" + inc, "-I" + C_DIR, "-o", tmp, os.path.join(C_DIR, name + ".c"), lib, xlib, "-lz", "-lm", "-lpthread"]
		p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=300)
		if p.returncode != 0:
			raise BuildError("failed to compile C driver %s:\n$ %s\n%s" % (name, " ".join(cmd), p.stdout.decode("utf-8", "replace")))
		os.replace(tmp, exe)
	return exe


def tree_products(root=None):
	"""Build products present in a source tree, as {relpath: mtime}."""
	import glob
	root = root or src_dir()
	out = {}
	for pat in TREE_PRODUCT_PATTERNS:
		for p in glob.glob(os.path.join(root, pat)):
			out[os.path.relpath(p, root)] = os.stat(p).st_mtime
	return out


if __name__ == "__main__":
	print(minibwa_bin())
