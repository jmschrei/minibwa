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
compiler changes. Nothing is written to the source tree.

| Variable | Default | Meaning |
|:---------|:--------|:--------|
| `MINIBWA_BIN` | build from source | binary under test |
| `MINIBWA_TEST_SRC_DIR` | repository root | source tree to build `minibwa` and `libminibwa.a` from |
| `CC` | `gcc` | compiler for the build and the C drivers |
| `MINIBWA_TEST_CACHE_DIR` | `$XDG_CACHE_HOME/minibwa-tests` or `~/.cache/minibwa-tests` | build cache |
| `MINIBWA_TEST_MAKE_JOBS` | 4 | `make -j` for the build |
| `MINIBWA_TEST_THREADS` | 1 | `-t` for tests that do not set it themselves |
| `MINIBWA_TEST_TIMEOUT` | 120 | seconds before one minibwa or driver process is killed |
| `MINIBWA_TEST_UPDATE_PINS` | unset | `1` rewrites pinned outputs (only on base minibwa) |
| `MBT_VERBOSE` | unset | C drivers also print passing checks |

`MINIBWA_BIN` replaces only the binary. The C drivers always link a
`libminibwa.a` built from `MINIBWA_TEST_SRC_DIR`, so to test a modified tree
point `MINIBWA_TEST_SRC_DIR` at it (and leave `MINIBWA_BIN` unset, or set it
to that tree's binary). The build cache is keyed by a hash of the Makefile,
the top-level `*.c`/`*.h` files, `mimalloc/` and the compiler's path and
version, and builds hold a file lock, so parallel workers and concurrent runs
are safe. Child processes run with `LC_ALL=C` and core dumps disabled.

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
`run_map(index, reads, args=(), threads=None, check=True)` and
`run_minibwa(argv)` return a `Result` with `.returncode .stdout .stderr
.stdout_bytes .sam() .paf() .check()`. `run_map` adds `-t
$MINIBWA_TEST_THREADS` unless `threads=` or a `-t N`/`-tN` argument is
given. Use `run_map(..., subcommand="mem")` for `minibwa mem`.

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
   the driver if needed (`$CC -std=gnu99 -O2 ... libminibwa.a -lz -lm
   -lpthread`), runs it and fails with its output unless every check passed.
   `H.c_cases(driver)` lists the cases.

Each failed check prints `not ok N - file:line: message`; a run ends with
`# mbt: CASE checks=N failures=M`. Pass inputs (index prefixes, FASTA files,
expected values) as arguments from Python; drivers should not build indexes.
`mbt_read_seqs()` reads a small plain FASTA/FASTQ file. Match minibwa's C
style (tabs, C99).

## Rules for expected values

* Derive expected values independently wherever practical: inputs small
  enough to work by hand, a short reference implementation inside the test
  (a naive suffix array or backward search, brute-force affine-gap
  alignment, NM/MD/cs recomputed from the genome), the simulator's recorded
  truth, the SAM specification. A test that only records the current output
  says so in its docstring ("pins the current output") and carries the
  `pinned` marker; keep these few.
* No test of behaviour base minibwa gets wrong. Such behaviour is reported
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

## Things to know about minibwa's output

* The `@PG` line carries `VN:` (`MB_VERSION`) and `CL:` with the full
  command line, including temporary paths. Compare or pin `strip_pg(text)`.
* stderr carries timings, CPU use, peak RSS and insert-size statistics.
  Never compare stderr wholesale; match the specific message.
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
* MD:Z: base minibwa omits the final `0` when the alignment ends with a
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
