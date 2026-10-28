# minibwa test suite

Unit and end-to-end tests of the mapping path: seeding and BWT lookups,
chaining, alignment, base-alignment tags, pairing and mate rescue, SAM/PAF
output, the `minibwa map` and `minibwa mem` command lines, and the
libminibwa API. Index building is not tested; tests build small indexes with
the binary under test. `TARGETS.md` lists every target and which test file
owns it.

## Running

Requirements: Python 3.8 or later, pytest, make, a C compiler and zlib.
samtools and pysam are optional; tests that use them are skipped when they
are missing.

```sh
python3 -m pytest tests                          # from the repository root
cd tests && python3 -m pytest                    # same
MINIBWA_TEST_THREADS=4 python3 -m pytest tests   # run minibwa with -t 4
python3 -m pytest tests -k pilot                 # a subset
```

The first run builds minibwa with the Makefile's own flags in a cache
directory (about 3 s); later runs reuse the build until a source file or the
compiler changes. Nothing is written to the source tree: minibwa processes
run in a temporary working directory.

### Compilers

`CC` selects the compiler for the build and for the C drivers, for example
`CC=gcc-11` or `CC=clang`. A clang installed in a conda environment (as
clang 18 was for the suite's own runs) does not find its zlib and OpenMP
runtime by itself; point the compiler, the linker and the loader at the
environment:

```sh
E=/path/to/conda/env    # holds bin/clang, include/zlib.h, lib/libz.so, lib/libomp.so
CC=$E/bin/clang CPATH=$E/include LIBRARY_PATH=$E/lib LD_LIBRARY_PATH=$E/lib \
	python3 -m pytest tests
```

`CPATH` and `LIBRARY_PATH` are read by the compiler for headers and
libraries (`zlib.h`, `-lz`, `-fopenmp`), and `LD_LIBRARY_PATH` lets the
built `minibwa` and the drivers load `libz.so` and `libomp.so` at run time.

### AddressSanitizer

`make asan=1` alone gives a binary that crashes at startup: mimalloc's
`malloc` override is called while the ASan runtime initializes. Build
without mimalloc instead, and without kalloc so that ASan sees every
allocation (an empty `MALLOC_O` keeps kalloc off; `mimalloc=0` would turn it
on), and compile the drivers with the sanitizer too:

```sh
MINIBWA_TEST_MAKE_ARGS="asan=1 MALLOC_O=" CC="gcc -fsanitize=address" \
	ASAN_OPTIONS=detect_leaks=0 python3 -m pytest tests
```

On minibwa r427 the whole suite passes this way with no ASan error. With
leak detection on, LeakSanitizer reports memory still allocated at exit
from `minibwa getref` and `minibwa fastmap`, and the tests that run them
fail on the exit code; it reports such memory too when `minibwa map`
returns early on an invalid `-R` line.

| Variable | Default | Meaning |
|:---------|:--------|:--------|
| `MINIBWA_BIN` | build from source | binary under test |
| `MINIBWA_TEST_SRC_DIR` | repository root | source tree to build `minibwa` and `libminibwa.a` from |
| `CC` | `gcc` | compiler for the build and the C drivers |
| `MINIBWA_TEST_CACHE_DIR` | `$XDG_CACHE_HOME/minibwa-tests` or `~/.cache/minibwa-tests` | build cache |
| `MINIBWA_TEST_MAKE_JOBS` | 4 | `make -j` for the build |
| `MINIBWA_TEST_MAKE_ARGS` | unset | extra `make` arguments for the build, e.g. `asan=1 MALLOC_O=` (part of the cache key) |
| `MINIBWA_TEST_THREADS` | 1 | `-t` for tests that do not set it themselves |
| `MINIBWA_TEST_TIMEOUT` | 120 | seconds before one minibwa or driver process is killed |
| `MINIBWA_TEST_UPDATE_PINS` | unset | `1` rewrites pinned outputs (only on current minibwa) |
| `MBT_VERBOSE` | unset | C drivers also print passing checks |

`MINIBWA_BIN` replaces only the binary. The C drivers always link a
`libminibwa.a` built from `MINIBWA_TEST_SRC_DIR`, so to test a modified tree
point `MINIBWA_TEST_SRC_DIR` at it (and leave `MINIBWA_BIN` unset, or set it
to that tree's binary). The build cache is keyed by a hash of the Makefile,
the top-level `*.c`/`*.h` files, `mimalloc/` and the compiler's path and
version, and builds hold a file lock, so parallel workers and concurrent runs
are safe. Child processes run with `LC_ALL=C`, in a fresh temporary working
directory unless a test passes `cwd=`, and with the core-file limit at 1
byte: no core is written, and on Linux a crash handler named in
`core_pattern` (such as apport, about a second per abort) is skipped, which
a limit of 0 would not do.

## Layout

| File | Contents |
|:-----|:---------|
| `minibwa_build.py` | locate or build the binary, `libminibwa.a` and C drivers out of tree |
| `minibwa_helpers.py` | genome and read simulation, FASTA/FASTQ I/O, running minibwa, SAM/PAF parsers, tag recomputation, pins, C driver runner |
| `conftest.py` | session fixtures |
| `c/minibwa_test.h` | check-and-report framework for C drivers |
| `c/drv_*.c` | C drivers, one binary each, compiled on demand |
| `pins/` | pinned outputs |
| `test_harness.py` | tests of the harness itself |
| `test_pilot.py` | whole-program sanity checks and one API case |
| `test_<area>_*.py` | tests per area of `TARGETS.md` (a, b, c, d) |

## Fixtures (conftest.py)

| Fixture | Scope | Value |
|:--------|:------|:------|
| `minibwa_bin` | session | path of the binary under test |
| `threads` | session | `MINIBWA_TEST_THREADS` |
| `index_dir` | session | a directory that lives for the session |
| `index_of` | session | `index_of(genome_or_fasta, meth=False)` -> index prefix, built once per genome |
| `chrm` | session | `.fasta`, `.genome` ({name: seq}), `.index`, `.read1`, `.read2` for the shipped `test/` data |
| `synth` | session | `.genome` (a `Genome`), `.fasta`, `.index` for the standard synthetic genome (`conftest.make_synth()`) |
| `tree_before` | session | build products in the source tree before the session |

## Helpers (minibwa_helpers.py)

Simulation, all deterministic for a given seed on every Python version
(the generator is splitmix64 in integer arithmetic; never use `random` or
`hash()` for inputs):

* `make_genome(seed, contig_lengths, repeat_families, repeat_len, repeat_copies, repeat_divergence, soft_mask_repeats, n_runs, n_run_len, low_complexity, low_complexity_len, short_contigs, short_contig_len)` -> `Genome` with `.names`, `[name]`, `.features` (`Feature(contig, start, end, kind, info)`, kind `repeat`/`N`/`lowcomp`), `.is_unique(contig, start, end, pad)`, `.write_fasta(path)`.
* `simulate_reads(genome, n, length, seed, sub_rate, indel_rate, max_indel, end_margin, n_rate, unique, contigs, strand, qual, name_prefix)` -> `[SimRead]`. A `SimRead` has `.name .seq .qual` (as sequenced) and the truth `.contig .pos .end .strand .cigar .edits .fwd_seq .nm`; indels stay `end_margin` bases from the read ends.
* `simulate_pairs(...)` -> `[SimPair]` (FR; `.r1 .r2 .contig .frag_start .frag_end .strand .insert`).
* `simulate_hic_pairs(..., inter_frac, min_distance, chimeric_frac)`: mates from independent loci; chimeric read 1 has two `.segments` (`Segment(contig, pos, end, strand, qstart, qend, cigar, edits)`).
* `simulate_long_reads(genome, n, length=(1000, 5000), ...)`.
* `Rng(seed)`: `u64 below randint random chance choice normal_int dna fork`.

I/O: `write_fastx/write_fasta/write_fastq(path, records)` (format and gzip
from the extension; gzip bytes are stable), `write_pairs(p1, p2, pairs,
suffix=False)`, `write_interleaved(path, pairs)`, `fastx_text`,
`read_fastx`, `read_fasta`. Records may be `SimRead`s or
`(name, seq[, qual[, comment]])` tuples.

Running: `index_genome(fa_or_genome, dir, meth=False)`;
`run_map(index, reads, args=(), threads=None, check=True, cwd=None)` and
`run_minibwa(argv, cwd=None)` return a `Result` with `.returncode .stdout
.stderr .stdout_bytes .sam() .paf() .check()`. `run_map` adds `-t
$MINIBWA_TEST_THREADS` unless `threads=` or a `-t N`/`-tN` argument is
given, and makes the index and read paths absolute. Both run the process
in a new temporary directory (removed afterwards) unless `cwd` is given,
so give other paths in `args` (such as `-o`) as absolute paths. Use
`run_map(..., subcommand="mem")` for `minibwa mem`. A process that aborts
returns `-6` without a core dump or crash-handler delay; no wrapper is
needed.

Parsing: `parse_sam(text, strict=True)` -> `SamFile` (`.hd .sq .rg .pg .co
.records .by_qname() .primaries()`); `SamRecord` fields `qname flag rname
pos mapq cigar rnext pnext tlen seq qual tags tag_types`, with flag
properties (`is_reverse`, `is_secondary`, `segment`, ...) and coordinates
(`pos0 ref_end0 left_clip unclipped_start0 query_length`). Tags are typed
(`i` int, `f` float, `A`/`Z`/`H` str, `B` `(subtype, [values])`); strict
parsing rejects duplicate tags, CIGAR/SEQ/QUAL length mismatches and
out-of-range fields. `parse_paf(text)` -> `[PafRecord]`.

Recomputation from the genome: `compute_nm`, `compute_md`, `compute_cs`
(`ref, pos0, seq, cigar`), `cigar_ref_span`, `cigar_query_span`,
`parse_cigar`, `revcomp`, `tlen_5prime`, `tlen_outer`, `placement_error(rec,
read)`. Generic checks returning a list of problems:
`sam_invariant_errors(sam, genome, reads)` and `pair_invariant_errors(sam)`.

Pins: `check_pinned(name, text)` compares with `pins/<name>` (gzip when the
name ends in `.gz`) and shows a diff; `strip_pg(sam)` drops the `@PG` line
and `sam_body(sam)` all header lines.

## Adding a C driver

1. Create `c/drv_<area><topic>.c` that includes `minibwa_test.h` and any of
   minibwa's headers (`minibwa.h`, or internal ones such as `mbpriv.h`,
   `bwt.h`, `ksw2.h`); the include path holds the source tree.
2. Write each case as `static int name(int argc, char *argv[])` using
   `MBT_CHECK(cond)`, `MBT_CHECKF(cond, fmt, ...)`, `MBT_CHECK_INT(a, op, b)`,
   `MBT_CHECK_STR(a, b)`, `MBT_REQUIRE(cond)` (stops the case) and
   `MBT_FAIL(fmt, ...)`. `argv` holds the arguments after the case name.
3. End the file with a `{name, function}` table terminated by `{0, 0}` and
   `MBT_MAIN(table)`.
4. From Python: `H.assert_c_case("drv_<...>", "case", arg1, ...)` compiles
   the driver if needed (`$CC -std=gnu99 -O2 ... libminibwa.a
   libmbt_extra.a -lz -lm -lpthread`), runs it and fails with its output
   unless every check passed. `H.c_cases(driver)` lists the cases.

`libminibwa.a` lacks some objects that `make minibwa` builds; the harness
archives `bseq.o` (`mb_bseq_*`) and `kthread.o` (`kt_for`, `kt_pipeline`)
from the same build into `libmbt_extra.a`, so a driver can call them
through `bseq.h` and `kthread.h` without including the `.c` file. Being an
archive, it adds only the members a driver needs, so a driver that does
include a `.c` file still links. `make` itself is not changed.

Each failed check prints `not ok N - file:line: message`; a run ends with
`# mbt: CASE checks=N failures=M`. Pass inputs (index prefixes, FASTA files,
expected values) as arguments from Python; drivers should not build indexes.
`mbt_read_seqs()` reads a small plain FASTA/FASTQ file. Match minibwa's C
style (tabs, C99).

Drivers that call only `minibwa.h` (`drv_pilot`, `drv_c_pe`) or a few
internal functions through their prototypes need at most call-site edits
when an internal signature changes. These read or build internal structs or
test internal components directly, and would need more than that if those
change: `drv_a_bwt` (`mb_bwt_t` `seq_len`/`primary`/`L2`, `l2b_t` and
`l2b_ctg_t` fields, `mb_sai_t`, `mb_anchor_t`, `mb_smem_entry_t`),
`drv_b_chain` and `drv_b_hits` (build `mb_anchor_t`, `mb_hit_t`,
`mb_extra_t` and `l2b_t` by hand), `drv_b_extra` (`mb_hit_t`/`mb_extra_t`,
the cs string after `cigar[n_cigar]`), `drv_b_ksw` (`ksw_extz_t`; vendored
ksw2), `drv_dkommon` (`l2b_t` fields, the `--meth` BWT text layout, and
`kt_for`/`kt_pipeline` semantics if kthread is replaced) and `drv_dapi`'s
`fmt` case (calls the internal `mb_format()` with an `mb_bseq1_t` and
`idx->l2b`). `drv_a_input` uses `mb_bseq1_t`, public `mb_opt_t` fields and
the public API.

## Rules for expected values

* Derive expected values independently wherever practical: inputs small
  enough to work by hand, a short reference implementation inside the test
  (a naive suffix array or backward search, brute-force affine-gap
  alignment, NM/MD/cs recomputed from the genome), the simulator's recorded
  truth, the SAM specification. A test that only records the current output
  says so in its docstring ("pins the current output") and carries the
  `pinned` marker; keep these few.
* No test of behaviour minibwa r427 gets wrong. Such behaviour is reported
  separately, not tested.
* No dependence on platform, compiler, Python version or thread count unless
  the test is explicitly about thread-count invariance (`thread_invariance`
  marker). minibwa builds with SSE on x86-64 and NEON on arm64, on Linux and
  macOS: do not depend on SIMD width, `ru_maxrss` units, `sed -i` dialects,
  `/proc`, or GNU-only tool options.
* Each test file runs fast; the whole suite should stay under about two
  minutes at `MINIBWA_TEST_THREADS=4` on a loaded machine. A minibwa process
  on a small input takes about 20 ms, most of it building the 10-mer cache at
  index load; prefer one run checked many ways over many runs.
* Do not pin implementation details that a speed optimization could change
  while keeping `minibwa map` output byte-identical: internal struct
  layouts, allocation counts, kalloc statistics, the order of internal calls,
  the order of intervals returned by `mb_seed_intv_batch`, the `hash`, `as`,
  `id` or `flt` fields of `mb_hit_t`. Pinning the bytes of `minibwa map`
  output on fixed inputs is fine and wanted.
* Never compare stderr wholesale. For a warning or an error, check that a
  message appears (a stderr line other than the `[M::...]` progress lines)
  and, where it carries one, that it names the offending option or value; do
  not match its wording: no minibwa document fixes it, and a faster
  implementation may reword it. Debug output (`--dbg-*`) may be parsed for
  its fields and numbers, but do not require a particular number of DP
  calls, anchors or trace lines.

## Things to know about minibwa's output

* The `@PG` line carries `VN:` (`MB_VERSION`) and `CL:` with the full
  command line, including temporary paths. Compare or pin `strip_pg(text)`.
* stderr carries timings, CPU use, peak RSS and insert-size statistics.
  Never compare stderr wholesale; match the specific message. The timing
  lines (`[M::...] mapped ...`, `Real time`) and the hidden `--dbg-*`
  traces are not asserted. The insert-size report of `mb_pestat` is
  undocumented, so `helpers_c.stderr_pestat` reads only its numbers (four
  candidate counts; quartiles, mean and standard deviation; bounds and
  percentages per orientation), not its wording.
* Output is byte-identical for every `-t` at a fixed `-K` (checked at
  `-t 1, 2, 4, 8` on the chrM data). For paired-end input whose insert size is
  inferred, the statistics are computed per batch, so `-K` changes the
  output; single-end output, and paired-end output with `-I`, do not depend
  on `-K`.
* A batch with fewer than 20 pairs uses the fixed insert-size statistics
  (`-I` or the defaults 400, 100, 50, 800) instead of inferring them. An
  orientation is inferred only from at least 20 unique pairs.
* Read names seed a hash that breaks ties between equally good hits.
  Renaming a read can move a read from a repeat to another copy; do not
  assert which copy a repetitive read lands on.
* A missing or unreadable index makes `kom_assert` call `abort()`: the
  Python return code is `-6` (SIGABRT) and stderr says
  `[E::main_map] failed to load the index. ABORT!`. A non-existent read file
  prints `ERROR: failed to open file ...` and writes only the header, but the
  exit status is 0; do not assert that exit status.
* Three adjacent reads with the same name (or two same-named pairs in
  two-file input) in paired-end mode fail an assertion and abort.
* MD:Z: minibwa r427 omits the final `0` when the alignment ends with a
  mismatch (`150T` instead of `150T0`). Compare with
  `md_from_minibwa(compute_md(...))`.
* TLEN is the 5'-to-5' template length (`tlen_5prime`); the SAM spec leaves
  the definition to the implementation. It differs from leftmost-to-rightmost
  (`tlen_outer`) for overlapping or dovetailed pairs.
* NM counts every non-ACGT base as a mismatch, as the SAM spec says. With
  `--meth`, NM does not count bisulfite conversions, but MD:Z and cs:Z do show
  them as mismatches.
* Supplementary records are hard-clipped (`-Y` soft-clips them). Secondary
  records have SEQ and QUAL `*` and soft clips, and are not written at all
  unless `--outn` is set; by default secondaries appear only in `XA:Z` (when
  at most `--xa` hits score above `--outs` times the best) and are counted in
  `n2:i`.
* PAF: names of paired reads get `/1` or `/2`; `tp:A:P` marks every primary
  chain, supplementary pieces included; unmapped reads get a line with `*`
  fields unless `-u`.
* A trailing `/1` or `/2` is removed from the names of paired reads in SAM.
* N in the reference is stored as a random base in the BWT and read back as
  N for alignment, so seeds can pass through reference N runs while the
  aligned N bases count as mismatches. N in a read breaks seeds.
* Reads shorter than `-l` (325) are aligned in short-read mode; longer reads
  in the default adaptive mode switch to long-read heuristics and a blended
  MAPQ. MAPQ is at most 60.
* `minibwa index` uses 4 threads unless `-t` is given; the helpers pass
  `-t 1` and set `OMP_NUM_THREADS`.
* The default build disables kalloc (`MB_F_NO_KALLOC`) unless built with
  `mimalloc=0`; `--kalloc=yes` turns it on.
* AS:i is the DP score (`mb_extra_t.dp_score`). It usually equals the score
  of the CIGAR, but keeps the cost of an indel that `mb_fix_cigar` slides to
  the read start and then drops (seen in (CTT)n tracts).
* `-O`/`-E` pairs get the documented cost min(O1+k*E1, O2+k*E2) only when
  O1 < O2 and E1 > E2 (reported in the issue linked from the pull request
  that added the tests).
* A size-1 SA interval at the sentinel's row is extended by one base that is
  not there: reads that run past the start or the end of the indexed text
  F + revcomp(F) can get SMEMs that do not exist (reported in the issue
  linked from the pull request that added the tests). Test reads there
  avoid the triggering base.
* Pair orientations: `mb_pestat` labels a same-strand pair FF or RR by read
  order, but pairing looks it up by coordinate order, so in an FF-only (or
  RR-only) library the pairs from reverse-strand fragments are never proper
  (reported in the issue linked from the pull request that added the
  tests). FR and RF libraries are not affected.

## Coverage

Line coverage of one run of the whole suite (450 tests, `MINIBWA_TEST_THREADS=4`)
on a copy of minibwa r427 built with `CC="gcc --coverage"` and `-O0` (gcc 13,
measured with gcovr 8.7). "CLI tests" counts what the `minibwa` processes
executed, "C drivers" what the drivers in `c/` executed (the two were kept
apart with `GCOV_PREFIX`), and "Both" their union. The mapping-path files
are minibwa's own sources reached by `minibwa map`, `minibwa mem`, the API
and the debugging subcommands; `bwt.c` and `l2bit.c` also hold
index-building code, which runs whenever a test builds an index. ksw2 and
mimalloc are vendored and listed separately; mimalloc is linked only into
the binary. `kseq.h` and `kalloc.h` hold only macros and declarations.

| File | CLI tests | C drivers | Both |
|:-----|----------:|----------:|-----:|
| `seed.c` | 90.1% (210/233) | 98.3% (229/233) | 99.1% (231/233) |
| `bwt.c` | 91.3% (409/448) | 82.1% (368/448) | 92.4% (414/448) |
| `l2bit.c` | 82.2% (277/337) | 49.0% (165/337) | 82.8% (279/337) |
| `lchain.c` | 99.3% (136/137) | 98.5% (135/137) | 99.3% (136/137) |
| `map-algo.c` | 84.6% (427/505) | 95.4% (482/505) | 99.4% (502/505) |
| `align.c` | 97.4% (668/686) | 75.1% (515/686) | 97.4% (668/686) |
| `pe.c` | 98.8% (422/427) | 71.4% (305/427) | 98.8% (422/427) |
| `cs.c` | 91.7% (111/121) | 100.0% (121/121) | 100.0% (121/121) |
| `format.c` | 96.1% (245/255) | 62.0% (158/255) | 96.5% (246/255) |
| `options.c` | 100.0% (90/90) | 100.0% (90/90) | 100.0% (90/90) |
| `map-main.c` | 99.8% (502/503) | 0.0% (0/503) | 99.8% (502/503) |
| `bseq.c` | 96.9% (94/97) | 94.8% (92/97) | 97.9% (95/97) |
| `kommon.c` | 84.4% (103/122) | 79.5% (97/122) | 90.2% (110/122) |
| `kthread.c` | 100.0% (68/68) | 100.0% (68/68) | 100.0% (68/68) |
| `kalloc.c` | 78.9% (90/114) | 78.9% (90/114) | 78.9% (90/114) |
| `main.c` | 99.2% (132/133) | 0.0% (0/133) | 99.2% (132/133) |
| `fastmap.c` | 100.0% (111/111) | 0.0% (0/111) | 100.0% (111/111) |
| `bwt.h` | 100.0% (6/6) | 100.0% (6/6) | 100.0% (6/6) |
| `l2bit.h` | 100.0% (4/4) | 100.0% (4/4) | 100.0% (4/4) |
| `mbpriv.h` | 65.8% (25/38) | 100.0% (38/38) | 100.0% (38/38) |
| `kommon.h` | 100.0% (5/5) | 100.0% (5/5) | 100.0% (5/5) |
| `bseq.h` | 77.8% (7/9) | 100.0% (9/9) | 100.0% (9/9) |
| `ksw2_extz2_sse.c` | 88.5% (154/174) | 100.0% (174/174) | 100.0% (174/174) |
| `ksw2_extd2_sse.c` | 94.2% (212/225) | 99.1% (223/225) | 100.0% (225/225) |
| `ksw2_ll_sse.c` | 94.0% (188/200) | 94.0% (188/200) | 94.0% (188/200) |
| `ksw2.h` | 90.1% (64/71) | 94.4% (67/71) | 94.4% (67/71) |
| **mapping-path .c (17 files)** | 93.3% (4095/4387) | 66.4% (2915/4387) | 96.1% (4217/4387) |
| **their headers** | 75.8% (47/62) | 100.0% (62/62) | 100.0% (62/62) |
| **ksw2 (vendored)** | 92.2% (618/670) | 97.3% (652/670) | 97.6% (654/670) |
| **mimalloc (vendored)** | 47.2% (2125/4501) | 0.0% (0/4501) | 47.2% (2125/4501) |

Lines no test reaches in the mapping path include: the `CG:B:I` branch of
`write_sam_cigar` (format.c:177-183, no caller sets it), the demotion of a
second primary overlapping the paired hit (pe.c:550-553), `mb_bseq_eof`
(no callers), closing already opened read files when a later one fails
(map-main.c:268), the failure path of `l2b_load_mmap`, and index-side code:
`mb_bwt_gen_sa` and `mb_bwt_load_raw` (bwt.c; the index is built through
libsais) and `l2b_save_pac`/`l2b_save_pac_meth` (l2bit.c).

Public names in `minibwa.h` that no test refers to: `MB_CIGAR_SOFTCLIP`,
`MB_CIGAR_HARDCLIP`, `MB_CIGAR_N_SKIP` and `MB_CIGAR_PADDING` (no hit's
CIGAR holds these operations; clips follow from `qs`/`qe`), the flag name
`MB_F_NO_KALLOC` (its effect is tested through `--kalloc` and
`mb_tbuf_init(no_kalloc)`), and the opaque struct tags `mb_idx_s` and
`mb_tbuf_s`. Every function declared there is called by a C driver.
