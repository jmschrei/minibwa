# Test targets

Every public name in `minibwa.h`, every function the mapping path reaches
(from `minibwa map`, `minibwa mem` and the libminibwa API), every
`minibwa map` option (`opt_str`, `long_options` and the usage text in
`map-main.c`), and the debugging subcommands that expose the seeding code.
Source revision: 74bdd48 (r427). Line numbers are for that revision.

Each row has an owner area. A writer ticks `[x]` in the last column when a
test covers the row, or writes `n/a` with a reason (for example "internal:
not observable without pinning implementation details", or "not reached").
Several rows can be covered by one test; name the test in Notes if useful.

## Area split

| Area | Scope | Test files it owns |
|:-----|:------|:-------------------|
| A | input and seeding: bseq.c, kseq parsing, batching (`-K`), l2bit.c lookups, bwt.c lookups, seed.c, `-k -c -K -y`, `--dbg-seed --dbg-anchor`, and the `fastmap`, `getref` and `bench` subcommands as oracles | `tests/test_a_*.py`, `tests/c/drv_a*.c`, `tests/pins/a_*` |
| B | chaining and alignment: lchain.c, hit processing and MAPQ in map-algo.c, align.c, ksw2, cs.c, options.c presets and the adaptive mode, `-A -B -O -E -s -w -W -g -m -p -N -l -x --long --adap --chain-only -b --eqx`, long reads | `tests/test_b_*.py`, `tests/c/drv_b*.c`, `tests/pins/b_*` |
| C | paired-end and Hi-C: pe.c, insert-size inference and `-I`, `-P`, `--rescue`, `--pe`, `-U`, `--hic`, `-5`, proper-pair flags and every mate field written by format.c, the paired branch of `mb_map_batch` | `tests/test_c_*.py`, `tests/c/drv_c*.c`, `tests/pins/c_*` |
| D | output, CLI and API: the rest of format.c (SAM/PAF records and header, `-R -H -u -Y --outn --outs --xa --xa-ratio -o -f -a`), `main()`/`main_map()`/`main_mem()` parsing and errors, `--version --help`, `--mmap`, `--meth`, `--kalloc`, thread and batch invariance (`-t`, `-K`), kthread, the libminibwa API and every name in minibwa.h | `tests/test_d_*.py`, `tests/c/drv_d*.c`, `tests/pins/d_*` |

Changes from the split in the brief, from reading the source:

* `--adap` is in B, not C. It toggles `MB_F_ADAP`, which only selects the
  read-length-dependent parameter scaling in `mb_opt_adap()` (options.c) and
  `mb_is_sr_mode()`; nothing in pe.c reads it. It belongs with `-x adap` and `-l`.
* map-algo.c is not in the brief's file list but holds the core of the mapping
  path (`mb_map_sai`, hit selection, MAPQ, the API). Hit processing and MAPQ
  go to B, the API functions to D, `mb_sync_high_cov` (only called by
  `mb_pair`) to C.
* The mate half of `mb_fmt_sam` (RNEXT, PNEXT, TLEN, flags 0x1/0x2/0x8/0x20/0x40/0x80,
  MC:Z, MQ:i, placement of an unmapped mate) goes to C, since the brief gives
  C "proper-pair flags, mate fields"; the rest of format.c is D.
* `--meth` touches seed.c, l2bit.c, align.c and pe.c. Its functions
  (`l2b_intv2cid_meth`, `l2b_meth_convert`, `l2b_meth_rev`,
  `mb_anchor_split_meth`) stay with D, which owns `--meth`, so one writer
  owns the BS-seq fixture.
* `--mmap` loaders (`mb_idx_load_mmap`, `l2b_load_mmap`,
  `mb_bwt_load_mmap`, `kom_mmap_file`) go to D with the option.
* `-K` appears twice in the brief. A owns how `-K NUM1[,NUM2]` cuts batches
  (bseq.c); D owns output invariance across `-t` and `-K`.
* `minibwa mem` reaches the same mapping path, so its options are rows (D).
* `fastmap`, `getref` and `bench` are debugging subcommands; they are A's
  windows onto `mb_bwt_smem`, `l2b_getseq` and the rank/SA functions.

## Entry points

| ID | Entry point | Source | Area | Notes | Done |
|:---|:------------|:-------|:-----|:------|:-----|
| E1 | `minibwa map` | map-main.c:435 `main_map` | D | options below; every area drives it | [ ] |
| E2 | `minibwa mem` | map-main.c:631 `main_mem` | D | bwa-mem option names mapped onto mb_opt_t; always writes MD:Z | [ ] |
| E3 | `mb_map()` (one read) | map-algo.c:670 | D | as api-test/ex-one.c; `mt` 0/1/2; NULL tbuf allowed | [ ] |
| E4 | `mb_map_batch()` | map-algo.c:695 | D (SE), C (PE) | as api-test/ex-batch.c; pairs reads 2i,2i+1 when `MB_F_PE` | [ ] |
| E5 | `minibwa fastmap` | fastmap.c:85 | A | prints SMEM intervals and positions; oracle for `mb_bwt_smem`/`mb_bwt_smem_batch` (`-b` > 1) | [ ] |
| E6 | `minibwa getref` | main.c:97 | A | prints the genome from .l2b via `l2b_getseq` (N for ambiguous) | [ ] |
| E7 | `minibwa bench` | main.c:133 | A | `-b 2a/sa/msa -n -v -p -1`; `-p` prints per-point rank/SA values | [ ] |
| E8 | `minibwa version`, `minibwa --help`, unknown command | main.c:55 | D | exit codes, stdout vs stderr | [ ] |

## Area A: input and seeding

| ID | Target | Source | Notes | Done |
|:---|:-------|:-------|:------|:-----|
| A01 | `mb_bseq_open` | bseq.c:43 | gz or plain via zlib; `-` reads stdin; NULL on failure | [ ] |
| A02 | `mb_bseq_close` | bseq.c:55 | | [ ] |
| A03 | `kseq2bseq` (static) | bseq.c:70 | U/u -> T/t; warning on empty name; qual kept unless PAF; comment kept with `-y` | [ ] |
| A04 | `kstrdup` (static) | bseq.c:62 | | [ ] |
| A05 | `mb_bseq_read` | bseq.c:85 | batch ends at `-K` NUM1 once >= 40000 reads, or at NUM2; in frag mode the mate of the last read joins the batch; parse-error warning | [ ] |
| A06 | `mb_bseq_read_frag` | bseq.c:133 | two files read in lockstep; unequal counts -> warning, extra records dropped | [ ] |
| A07 | `mb_bseq_eof` | bseq.c:163 | not reached (no callers) | [ ] |
| A08 | `mb_qname_len`, `mb_qname_same` | bseq.h:26 | `/<digit>` suffix ignored when matching mates | [ ] |
| A09 | kseq parser (`KSEQ_INIT2` in bseq.c, kseq.h) | kseq.h | multi-line FASTA, FASTQ, mixed FASTA/FASTQ, comments, empty sequence | [ ] |
| A10 | read grouping in `worker_pipeline` step 0 | map-main.c:112 | adjacent same-name reads form a pair (interleaved); trailing `/1` `/2` trimmed from paired names; `assert(i - j <= 2)` | [ ] |
| A11 | sub-batching in `worker_pipeline` (`sb_len`, `sb_seq`) | map-main.c:162 | output must not depend on where sub-batches split | [ ] |
| A12 | `worker_for_se_batch` | map-main.c:36 | nt4 encoding, `mb_seed_intv_batch`, `mb_opt_adap` per read | [ ] |
| A13 | `mb_seed_intv` | seed.c:24 | two passes: SMEMs, then sub-SMEMs (min length half the SMEM, occurrence > parent) for SMEMs >= 2k with occ <= max_sub_occ | [ ] |
| A14 | `mb_seed_intv_batch` | seed.c:56 | same intervals as `mb_seed_intv`, order may differ | [ ] |
| A15 | `mb_seed_sort_dedup` (static) | seed.c:109 | | [ ] |
| A16 | `mb_anchor` | seed.c:268 | SA lookup; intervals with occ > `-c` sampled to `-c` positions; returns seed ratio (1.0 unless meth) | [ ] |
| A17 | `process_batch` (static) | seed.c:176 | SA position -> contig/strand via `l2b_intv2cid`; anchors spanning contigs dropped | [ ] |
| A18 | `mb_anchor_dedup` (static) | seed.c:143 | drops anchors contained in another on the same diagonal | [ ] |
| A19 | `mb_anchor_sort` | seed.c:341 | sort by concatenated target position (also used by B's re-chaining) | [ ] |
| A20 | radix sorts `mb_sai0`, `mb_sais`, `mb_saii`, `mb_anchor` | seed.c:9-18 | | [ ] |
| A21 | `mb_bwt_load`, `read_huge` | bwt.c:644, 593 | loads .mbw; failure -> NULL -> `kom_assert` abort in map | [ ] |
| A22 | `mb_bwt_destroy` | bwt.c:34 | | [ ] |
| A23 | `mb_bwt_cache`, `mb_bwt_count_kmer` | bwt.c:581, 548 | 10-mer interval table built at load; each entry equals the interval found by backward search | [ ] |
| A24 | `mb_bwt_set_intv` (inline) | bwt.h:75 | | [ ] |
| A25 | `mb_bwt_rank2a`, `mb_bwt_rank1a`, `mb_bwt_rank11`, `seek_block`, `rank_aux1`, `rank_aux4`, `mb_bwt_block_prefetch` | bwt.c:121-235 | occurrence counts equal a naive count over the BWT string | [ ] |
| A26 | `mb_bwt_extend` | bwt.c:237 | bidirectional extension; intervals equal naive suffix-array intervals of both strands | [ ] |
| A27 | `mb_bwt_back` (static), `mb_bwt_smem` | bwt.c:253, 280 | SMEMs equal a brute-force SMEM search; N breaks seeds; `min_occ` | [ ] |
| A28 | `mb_bwt_smem_batch`, `tq_init`, `tq_push`, `tq_shift`, `se_one_step_back` | bwt.c:316-460 | same intervals as `mb_bwt_smem` for any batch size | [ ] |
| A29 | `bwt_invPsi`, `mb_bwt_sa`, `mb_bwt_sa_batch` | bwt.c:465-546 | SA values equal a naive suffix array; batch equals one-by-one | [ ] |
| A30 | index-side bwt.c: `mb_bwt_init`, `mb_bwt_init_from_raw`, `mb_bwt_data_len`, `bwt_gen_cnt_table`, `mb_bwt_gen_sa`, `mb_bwt_save`, `mb_bwt_load_raw` | bwt.c | out of scope (index building); fixture only | n/a |
| A31 | `l2b_load` | l2bit.c:275 | contig names, lengths, offsets, ambiguous intervals | [ ] |
| A32 | `l2b_destroy` | l2bit.c:231 | | [ ] |
| A33 | `l2b_pos2cid` (static), `l2b_intv2cid` | l2bit.c:9, 25 | -1 spans two contigs, -2 spans the strand boundary, -3 out of range; reverse-strand offsets | [ ] |
| A34 | `l2b_getseq` | l2bit.c:52 | codes 0-4, 4 for reference ambiguity codes (N and IUPAC); clamps to the contig | [ ] |
| A35 | `l2b_getambi` | l2bit.c:91 | | [ ] |
| A36 | `l2b_get0`, `l2b_seq_prefetch` (inline) | l2bit.h:52-61 | | [ ] |
| A37 | index-side l2bit.c: `l2b_import`, `l2b_format_seq`, `l2b_add_seq`, `l2b_collate_str`, `l2b_save`, `l2b_save_pac`, `l2b_save_pac_meth`, `l2b_c2t`, `l2b_g2a` | l2bit.c | out of scope; note N in the FASTA becomes a random base (seed 11) in the BWT and 4 in `l2b_getseq`; lowercase is a soft mask only | n/a |
| A38 | `minibwa fastmap` internals: `main_fastmap`, `batch_smem`, `write_intv`, `process_batch` | fastmap.c | `-l -s -w -b` | [ ] |
| A39 | `main_getref` | main.c:97 | | [ ] |
| A40 | `main_bench` | main.c:133 | | [ ] |
| A41 | `kom_nt4_table` | kommon.c:162 | A/C/G/T any case -> 0-3, everything else -> 4 (U is converted to T before) | [ ] |

### Area A options

| ID | Option | Usage text | Notes | Done |
|:---|:-------|:-----------|:------|:-----|
| A50 | `-k INT` | min seed length [19] | `mb_opt_t.min_len`; also the sub-SMEM threshold and `chn_pen_gap` scale | [ ] |
| A51 | `-c NUM` | max seed occurrences [250] | `max_occ`; `kom_parse_num` suffixes; sampling of repetitive seeds | [ ] |
| A52 | `-K NUM1[,NUM2]` | process NUM1-NUM2 bp of query sequences in a batch [100m,1g] | `mb_size`, `max_mb_size`; NUM2 < NUM1 is raised to NUM1 | [ ] |
| A53 | `-y` | copy FASTA/Q comments to output | `MB_F_COPY_COMMENT`; appended verbatim to SAM and PAF lines | [ ] |
| A54 | input forms | `<in.fastq>` [read2] | plain/gzip, FASTA/FASTQ, two files, interleaved, stdin `-`, mixed SE and PE | [ ] |
| A55 | `--dbg-seed` | (hidden) | `SD` lines on stderr | [ ] |
| A56 | `--dbg-anchor` | (hidden) | `AC` lines on stderr | [ ] |

### Area A mb_opt_t fields

| ID | Field | Default | CLI | Done |
|:---|:------|:--------|:----|:-----|
| A60 | `min_len` | 19 | `-k` | [ ] |
| A61 | `max_sub_occ` | 10 | none | [ ] |
| A62 | `max_occ` | 250 | `-c` | [ ] |
| A63 | `sb_len` | 1000000 | none | [ ] |
| A64 | `sb_seq` | 24 | none | [ ] |
| A65 | `mb_size` | 100000000 (sr/adap), 500000000 (lr) | `-K` | [ ] |
| A66 | `max_mb_size` | 1000000000 | `-K` | [ ] |

## Area B: chaining and alignment

| ID | Target | Source | Notes | Done |
|:---|:-------|:-------|:------|:-----|
| B01 | `mb_lchain_dp` | lchain.c:165 | chaining DP; `max_skip`, `max_iter`, `bw`, `max_gap`, `min_sc` | [ ] |
| B02 | `comput_sc` (static) | lchain.c:133 | anchor-pair score with linear+log gap cost | [ ] |
| B03 | `mb_chain_backtrack` (static) | lchain.c:43 | | [ ] |
| B04 | `mb_chain_bk_end` (static) | lchain.c:25 | | [ ] |
| B05 | `compact_a` (static) | lchain.c:94 | chains reordered by target position | [ ] |
| B06 | `mb_map_sai` | map-algo.c:578 | seeds -> anchors -> chains -> re-chaining with `-W` for long reads -> alignment -> MAPQ | [ ] |
| B07 | `mb_cal_mblen`, `mb_cal_fuzzy_len` (static) | map-algo.c:133, 151 | chain-level mlen/blen (PAF columns 10-11 with `--chain-only`) | [ ] |
| B08 | `mb_hit_set_coor` (static) | map-algo.c:156 | | [ ] |
| B09 | `mb_cal_high_cov` | map-algo.c:175 | query fraction covered by seeds with occ > `-c` | [ ] |
| B10 | `mb_gen_hit` | map-algo.c:210 | ties broken by a hash of the read name (`mb_hash_str`) | [ ] |
| B11 | `mb_split_hit` | map-algo.c:247 | Z-drop split -> supplementary | [ ] |
| B12 | `mb_sync_hits` | map-algo.c:266 | | [ ] |
| B13 | `update_sub` (static), `mb_set_parent` | map-algo.c:294, 315 | primary/secondary by query overlap (`mask_level`) | [ ] |
| B14 | `mb_set_sam_pri` | map-algo.c:366 | first primary, or the smallest query start with `-5` (C owns `-5`) | [ ] |
| B15 | `mb_select_sub` | map-algo.c:383 | `-p` ratio and `-N` cap on secondaries | [ ] |
| B16 | `mb_hit_sort` | map-algo.c:407 | by DP score | [ ] |
| B17 | `mb_filter_hits` | map-algo.c:435 | `-m` on mlen and `-s` on dp_max | [ ] |
| B18 | `mb_squeeze_a` | map-algo.c:454 | | [ ] |
| B19 | `mb_set_inv_mapq` (static), `mb_set_mapq` | map-algo.c:478, 504 | short-read (bwa-mem) and long-read (minimap2) formulas, blended in adaptive mode; 0-60 | [ ] |
| B20 | radix sorts `mb128x`, `mb64` | map-algo.c:10-13 | | [ ] |
| B21 | `update_max_zdrop` (static), `mm_test_zdrop` (static) | align.c:11, 26 | | [ ] |
| B22 | `mb_fix_cigar` (static) | align.c:70 | removes leading/adjacent I/D artefacts | [ ] |
| B23 | `mm_update_cigar_eqx` (static) | align.c:148 | `--eqx`: `=`/`X`; N counts as `X` | [ ] |
| B24 | `mb_update_extra` | align.c:219 | mlen/blen/n_ambi, dp_max; writes cs/ds/MD | [ ] |
| B25 | `mb_enlarge_cigar` (static), `mb_append_cigar` | align.c:285, 300 | merges same ops at the junction | [ ] |
| B26 | `mb_min_int32`, `max_bw_from_mm` (static) | align.c:316, 321 | | [ ] |
| B27 | `mb_align_pair` (static) | align.c:329 | ungapped shortcuts, bandwidth cap, `max_sw_mat` skip, affine vs dual-affine choice | [ ] |
| B28 | `collect_long_gaps`, `mm_filter_bad_seeds`, `mm_filter_bad_seeds_alt`, `mm_fix_bad_ends`, `mb_max_stretch` (static) | align.c:392-546 | long-read anchor cleanup; short-read longest stretch | [ ] |
| B29 | `mb_align1` (static) | align.c:548 | left/right extension with end bonus, gap filling, Z-drop | [ ] |
| B30 | `mb_align1_inv` (static) | align.c:759 | inversion between two split hits (long reads) | [ ] |
| B31 | `mb_insert_reg`, `mb_count_gaps`, `mb_event_identity`, `mb_recal_max_dp` (static), `mb_update_dp_max` | align.c:817-892 | long-read rescoring | [ ] |
| B32 | `mb_align_skeleton` | align.c:894 | | [ ] |
| B33 | `ksw_extz2_sse` | ksw2_extz2_sse.c:16 | affine-gap extension/global alignment; compare with a brute-force affine DP | [ ] |
| B34 | `ksw_extd2_sse`, `ksw_i8x4_to_i32x4` | ksw2_extd2_sse.c:16, 27 | dual affine gap (`-O a,b -E c,d`) | [ ] |
| B35 | `ksw_ll_qinit`, `ksw_ll_u8_core`, `ksw_ll_i16_core`, `ksw_ll_i16`, `ksw_le_u8`, `ksw_max_u8`, `ksw_le_epi16`, `ksw_max_i16` | ksw2_ll_sse.c | local SW used by mate rescue (C) and inversions | [ ] |
| B36 | `ksw_push_cigar`, `ksw_backtrack`, `ksw_reset_extz`, `ksw_apply_zdrop`, `ksw_gen_nt4_mat` (inline) | ksw2.h:124-219 | scoring matrix: `-A`, `-B`, `b_ts`, `b_ambi`, meth columns | [ ] |
| B37 | `ksw_extz`, `ksw_extd`, `ksw_exts2_sse` | ksw2.h | declared, not compiled | n/a |
| B38 | `alloc_tmp`, `write_indel_ds` (static) | cs.c:7, 18 | | [ ] |
| B39 | `mb_write_cs_ds` | cs.c:47 | `-b cs` short cs; `-b ds` with `[..]` repeat context around indels | [ ] |
| B40 | `mb_write_MD` | cs.c:129 | `-b MD`; base writes no final `0` after a terminal mismatch (SAMtags regex wants one) | [ ] |
| B41 | `mb_opt_reset` (static), `mb_opt_preset` | options.c:5, 54 | sr/adap/lr values; unknown preset -> -1 | [ ] |
| B42 | `mb_opt_adap` | options.c:89 | `x1-(x1-x0)*2^(-(max(L,100)-100)/1900)` for -g -w -W -N -m -s (man page) | [ ] |
| B43 | `mb_log2`, `mb_hash64`, `mb_hash_str`, `mb_seq_rev`, `mb_is_sr_mode` (inline) | mbpriv.h:105-142 | output depends on read names through `mb_hash_str` | [ ] |
| B44 | `mb_dbg_seed`, `mb_dbg_anchor` (static) | map-algo.c:556, 565 | debug output (rows A55, A56) | [ ] |

### Area B options

| ID | Option | Usage text | Notes | Done |
|:---|:-------|:-----------|:------|:-----|
| B50 | `-A INT` | matching score [2] | `a` | [ ] |
| B51 | `-B INT` | mismatching openalty [8] | `b` (sic: "openalty" in usage) | [ ] |
| B52 | `-O INT1[,INT2]` | gap open penalty [12,23] | one value sets both | [ ] |
| B53 | `-E INT1[,INT2]` | gap extension penalty [2,1] | one value sets both; equal pairs select `ksw_extz2_sse` | [ ] |
| B54 | `-s INT` | suppress alignment with DP score lower than INT*{-A} [30] | `min_dp_max`; adaptive | [ ] |
| B55 | `-w NUM` | bandwidth [100] | `bw`; adaptive | [ ] |
| B56 | `-W NUM` | long bandwidth (for long reads or the adaptive mode) [20000] | `bw_long` | [ ] |
| B57 | `-g NUM` | max gap size, controlling extension and chain breaking [100] | `max_gap`; adaptive | [ ] |
| B58 | `-m INT` | min chaining score [25] | `min_chain_score`; adaptive | [ ] |
| B59 | `-p FLOAT` | min secondary-to-primary score ratio [0.5] | `pri_ratio` | [ ] |
| B60 | `-N INT` | retain at most INT secondary alignments [50] | `best_n`; adaptive | [ ] |
| B61 | `-l NUM` | treat reads <NUM as short reads in the default adaptive mode [325] | `max_sr_len` | [ ] |
| B62 | `-x STR` | preset (sr, lr or adap for mixed short/long reads) [adap] | applied before all other options regardless of position | [ ] |
| B63 | `--long[=yes/no]` | (man page) | `MB_F_LONG` | [ ] |
| B64 | `--adap=yes/no` | (man page) | `MB_F_ADAP`; invalid value -> warning | [ ] |
| B65 | `--chain-only` | perform chaining only without base alignment; force -fP | `MB_F_NO_ALN` | [ ] |
| B66 | `-b STR` | output a base alignment tag: cs, ds or MD [] | `md` accepted; invalid value -> warning and cs | [ ] |
| B67 | `--eqx` | (man page lists no text; README none) | `MB_F_EQX` | [ ] |
| B68 | `--dbg-aln-seq` | (hidden) | DP inputs on stderr | [ ] |
| B69 | `--dbg-an-pos` | (hidden) | `AF`/`AD` lines on stderr | [ ] |

### Area B mb_opt_t fields

| ID | Field | Default (adap/sr; lr) | CLI | Done |
|:---|:------|:----------------------|:----|:-----|
| B70 | `bw`, `bw_long` | 100, 20000; 500, 20000 | `-w`, `-W` | [ ] |
| B71 | `max_gap` | 100; 5000 | `-g` | [ ] |
| B72 | `max_sr_len` | 325 | `-l` | [ ] |
| B73 | `max_chain_skip`, `max_chain_iter` | 25, 5000 | none | [ ] |
| B74 | `min_chain_score` | 25; 40 | `-m` | [ ] |
| B75 | `chain_gap_scale` | 0.8 | none | [ ] |
| B76 | `mask_level`, `mask_len` | 0.5, INT_MAX | `mem -D` | [ ] |
| B77 | `pri_ratio` | 0.5 | `-p` | [ ] |
| B78 | `best_n` | 50; 5 | `-N` | [ ] |
| B79 | `a`, `b`, `b_ts`, `b_ambi` | 2, 8, 0, 1 | `-A`, `-B` | [ ] |
| B80 | `q`, `q2`, `e`, `e2` | 12, 23, 2, 1 | `-O`, `-E` | [ ] |
| B81 | `end_bonus` | 10; -1 | none | [ ] |
| B82 | `min_dp_max` | 30; 50 | `-s` | [ ] |
| B83 | `zdrop`, `zdrop_inv` | 80, 80; 400, 240 | none | [ ] |
| B84 | `min_ksw_len` | 20; 200 | none | [ ] |
| B85 | `max_sw_mat` | 100000000 | none | [ ] |
| B86 | `MB_F_LONG`, `MB_F_ADAP`, `MB_F_EQX`, `MB_F_NO_ALN`, `MB_F_WRITE_CS`, `MB_F_WRITE_DS`, `MB_F_WRITE_MD` | | flags above | [ ] |

## Area C: paired-end and Hi-C

| ID | Target | Source | Notes | Done |
|:---|:-------|:-------|:------|:-----|
| C01 | `mb_insert_dir` (static) | pe.c:10 | orientation code FF/FR/RF/RR and 5'-to-5' distance | [ ] |
| C02 | `mb_pair_score` (static) | pe.c:19 | | [ ] |
| C03 | `mb_select_unique_se` (static) | pe.c:32 | one primary with MAPQ >= 10 | [ ] |
| C04 | `mb_hit_sum_score` (static) | pe.c:41 | | [ ] |
| C05 | `mb_pestat` | pe.c:68 | per batch, >= 20 pairs, quartile bounds; stderr report | [ ] |
| C06 | `mb_pair_hits` (static) | pe.c:145 | | [ ] |
| C07 | `mb_ungap` (static) | pe.c:214 | 7-mer Hough filter before rescue | [ ] |
| C08 | `mb_matesw_align` (static) | pe.c:254 | | [ ] |
| C09 | `mb_matesw_core` (static), `mb_matesw` (static) | pe.c:324, 385 | rescue window from the insert bounds; `--rescue` candidates; `-U` | [ ] |
| C10 | `mb_pair` | pe.c:466 | pairing, pair MAPQ, primary lifting, chimeric hits, `-5` re-selection | [ ] |
| C11 | `worker_for_pe`, PE branch of `worker_pipeline` step 1 | map-main.c:96, 177 | fewer than 20 pairs in a batch or `-I` -> fixed stats (avg, std, lo, hi); else `mb_pestat` | [ ] |
| C12 | `set_ins_size` | map-main.c:424 | `-I` defaults: std = round(avg/10), hi = avg+4std, lo = max(avg-4std, 1) | [ ] |
| C13 | `mb_sync_high_cov` | map-algo.c:201 | only from `mb_pair` | [ ] |
| C14 | `mb_set_sam_pri` with `is_primary5` | map-algo.c:366 | `-5` | [ ] |
| C15 | `mb_fmt_sam` mate section | format.c:205-268, 296-302 | 0x1/0x2/0x8/0x20/0x40/0x80; RNEXT/PNEXT; TLEN; unmapped mate placed at the mapped mate; MC:Z and MQ:i | [ ] |
| C16 | `get_sam_pri` | format.c:158 | mate primary lookup | [ ] |
| C17 | paired branch of `mb_map_batch` | map-algo.c:754 | reads 2i and 2i+1 with the fixed (`-I`) stats; meth G2A on odd reads | [ ] |
| C18 | PAF `/1` `/2` name suffix | format.c:22 | | [ ] |
| C19 | `mb_hit_t.proper_pair`, `mb_hit_t.rescued` | minibwa.h | through `mb_map_batch` with `MB_F_PE` | [ ] |

### Area C options

| ID | Option | Usage text | Notes | Done |
|:---|:-------|:-----------|:------|:-----|
| C50 | `-P` | skip pairing and mate rescue | `MB_F_NO_PAIRING` | [ ] |
| C51 | `--rescue=INT` | mate rescue for up to INT candidates; 0 to skip rescue [10] | `max_rescue` | [ ] |
| C52 | `-I INT[,INT[,INT[,INT]]]` | mean, stddev, max and min of isize distribution [inferred] | `MB_F_PE_PREDEF` | [ ] |
| C53 | `--pe=yes/no` | (man page) | `MB_F_PE`; lr preset disables | [ ] |
| C54 | `-U INT` | (hidden) | `pen_unpair` [17] | [ ] |
| C55 | `--hic` | map Hi-C reads; equivalent to option -5P | | [ ] |
| C56 | `-5` | take the alignment with the smallest query position as primary | `MB_F_PRIMARY5` | [ ] |
| C57 | `--dbg-aln-pe` | (hidden) | rescue DP on stderr | [ ] |
| C58 | `--dbg-qname` | (hidden) | `QN`/`QP` lines on stderr | [ ] |

### Area C mb_opt_t fields

| ID | Field | Default | CLI | Done |
|:---|:------|:--------|:----|:-----|
| C60 | `max_pe_ins` | 10000 | none | [ ] |
| C61 | `max_rescue` | 10 | `--rescue` | [ ] |
| C62 | `pen_unpair` | 17 | `-U` | [ ] |
| C63 | `pe_avg`, `pe_std`, `pe_lo`, `pe_hi` | 400, 100, 50, 800 | `-I` | [ ] |
| C64 | `MB_F_PE`, `MB_F_PE_PREDEF`, `MB_F_PRIMARY5`, `MB_F_NO_PAIRING` | | flags above | [ ] |

## Area D: output, CLI and API

| ID | Target | Source | Notes | Done |
|:---|:-------|:-------|:------|:-----|
| D01 | `write_tags` (static) | format.c:13 | NM = blen - mlen + n_ambi; AS, ms, md | [ ] |
| D02 | `mb_fmt_paf` | format.c:19 | 12 columns; tp, s1, cm, s2 (primary only), NM, AS, ms, md, cg, cs/ds/MD, comment; unmapped line | [ ] |
| D03 | `mb_escape` | format.c:52 | `\t \n \r \\` | [ ] |
| D04 | `sam_write_rg_line` (static) | format.c:69 | `-R` must start with @RG, have ID:, no literal tab, ID < 256 chars | [ ] |
| D05 | `mb_fmt_sam_hdr` | format.c:105 | @HD VN:1.6 SO:unsorted GO:query; @SQ per contig; @RG; @PG with VN and CL | [ ] |
| D06 | `str_enlarge`, `str_copy`, `sam_write_sq` (static) | format.c:129-156 | SEQ reverse-complemented and QUAL reversed on the reverse strand | [ ] |
| D07 | `write_sam_cigar` (static) | format.c:168 | soft clips; hard clips on supplementary unless `-Y`; the `CG:B:I` branch is not reached | [ ] |
| D08 | `mb_fmt_sam` (non-mate parts) | format.c:196 | flags 0x4/0x10/0x100/0x800; SEQ `*` for secondaries; SA:Z; XA:Z and n2:i (`--xa`, `--outs`); RG:Z; comment | [ ] |
| D09 | `mb_format` | format.c:355 | SAM vs PAF | [ ] |
| D10 | output step of `worker_pipeline` | map-main.c:194 | `--outn` and `--outs` filtering of secondaries; `-u`; input order kept | [ ] |
| D11 | `mb_open_bseqs`, `mb_map_file` | map-main.c:258, 276 | `-o`; an unreadable read file prints an error but `main_map` still exits 0 | [ ] |
| D12 | `mb_insert_hdr` (static) | map-main.c:305 | `-H @...` (escapes) or `-H file` (lines `@XX\t...`); only the last `-H` is kept | [ ] |
| D13 | `usage_map` | map-main.c:360 | `--help` -> stdout, exit 0; missing arguments -> stderr, exit 1 | [ ] |
| D14 | `yes_or_no` (static) | map-main.c:411 | yes/y/no/n; other values warn | [ ] |
| D15 | `main_map` | map-main.c:435 | two-pass parsing (preset first); unknown option / missing argument / unknown preset -> exit 1; index load failure -> abort | [ ] |
| D16 | `usage_mem`, `main_mem` | map-main.c:578, 631 | see D70-D73 | [ ] |
| D17 | `main` dispatch | main.c:55 | `version`, `--help`, no arguments, unknown command; `[M::main]` stderr lines | [ ] |
| D18 | `kom_parse_num` | kommon.c:21 | k/m/g suffixes, decimals (`1.5k`) | [ ] |
| D19 | `kom_sprintf_lite`, `km_sprintf_lite`, `kom_sprintf_lite_core`, `str_enlarge`, `str_copy` | kommon.c:47-155 | every SAM/PAF byte goes through these | [ ] |
| D20 | `kom_strdup`, `kom_panic` | kommon.c:11, 33 | panic -> abort (SIGABRT) | [ ] |
| D21 | `kt_for`, `kt_pipeline`, `steal_work`, `ktf_worker`, `ktp_worker` | kthread.c | output identical for any `-t` | [ ] |
| D22 | kalloc.c (`km_init`, `kmalloc`, ...) | kalloc.c | `--kalloc=yes` output identical; internals not tested | [ ] |
| D23 | `mb_idx_load` | map-algo.c:19 | NULL for a missing prefix; `is_meth` picks .meth.mbw | [ ] |
| D24 | `mb_idx_load_mmap`, `l2b_load_mmap`, `mb_bwt_load_mmap`, `kom_mmap_file`, `kom_munmap` | map-algo.c:44, l2bit.c:332, bwt.c:681, kommon.c:348 | `--mmap`, `--mmap=lite`: same output as without | [ ] |
| D25 | `mb_idx_destroy` | map-algo.c:69 | NULL-safe | [ ] |
| D26 | `mb_idx_ctg_name`, `mb_idx_ctg_len` | map-algo.c:77, 82 | NULL / -1 out of range | [ ] |
| D27 | `mb_tbuf_init`, `mb_tbuf_km`, `mb_tbuf_destroy`, `mb_tbuf_reset` | map-algo.c:95-130 | results independent of the buffer and of kalloc | [ ] |
| D28 | `mb_map` | map-algo.c:670 | same hits as `minibwa map` for the read; `mt` != 0 on a non-meth index -> 0 hits; ASCII or 0-3 encoded input | [ ] |
| D29 | `mb_map_batch` (single-end) | map-algo.c:695 | same hits as `mb_map`; `n_seq` <= 0 -> NULL; sub-batching by `sb_len`/`sb_seq` | [ ] |
| D30 | `l2b_intv2cid_meth` | l2bit.c:36 | four converted copies | [ ] |
| D31 | `l2b_meth_convert`, `l2b_meth_rev` | l2bit.c:79, l2bit.h:63 | | [ ] |
| D32 | `mb_anchor_split_meth` (static) | seed.c:212 | splits anchors at non-conversion mismatches; seed ratio | [ ] |
| D33 | meth paths in `worker_for_se_batch`, `mb_map`, `mb_map_batch`, `mb_matesw`, `ksw_gen_nt4_mat` | | read 1 C->T, read 2 G->A; conversions score as matches and do not count in NM | [ ] |

### Area D options

| ID | Option | Usage text | Notes | Done |
|:---|:-------|:-----------|:------|:-----|
| D50 | `-f` | output PAF (SAM by default) | `MB_F_PAF` | [ ] |
| D51 | `-a` | (hidden) | clears `MB_F_PAF` (SAM) | [ ] |
| D52 | `-t INT` | number of worker threads [1] | | [ ] |
| D53 | `-R STR` | SAM read group line in a format like '@RG\tID:foo\tSM:bar' [] | error paths exit 1 | [ ] |
| D54 | `-o FILE` | output file name [stdout] | `-` is stdout | [ ] |
| D55 | `-u` | don't output unmapped reads | `MB_F_NO_UNMAP` | [ ] |
| D56 | `--outn=NUM` | output up to {NUM,-N} secondary alignments [0] | | [ ] |
| D57 | `--outs=FLOAT` | output a secondary hit if score at least FLOAT*bestScore [0.8] | also limits XA | [ ] |
| D58 | `--xa-ratio=FLOAT` | (hidden) | same as `--outs` | [ ] |
| D59 | `--xa=NUM` | if <=NUM hits with score >80% of the best hit, output them to XA [5] | 0 disables XA and n2 | [ ] |
| D60 | `-Y` | use soft clipping for supplementary alignments | `MB_F_SUPP_SOFT` | [ ] |
| D61 | `-H STR` | if STR starts with @, insert to header; or insert lines in file STR [] | | [ ] |
| D62 | `--mmap[=lite]` | load the index via memory mapped files (slower mapping) [] | | [ ] |
| D63 | `--meth` | map *directional* bisulfite sequencing reads | needs `index --meth`; else abort | [ ] |
| D64 | `--kalloc=yes/no` | (hidden) | `MB_F_NO_KALLOC` | [ ] |
| D65 | `--version` | print version number | `MB_VERSION` on stdout, exit 0 | [ ] |
| D66 | `--help` | print this help message | | [ ] |
| D67 | `-t`/`-K` invariance | | byte-identical output for `-t 1..N`; for `-K` with SE input or with `-I` (PE stats are per batch otherwise) | [ ] |
| D68 | errors | | unknown option, missing argument, unknown preset, missing index (abort, -6), unreadable reads, `-R` errors | [ ] |
| D70 | `mem -t -k -w -c -D -W -m -S -P` | usage_mem "Algorithm" | `-W` is min chaining score, `-m` is rescue count, `-D` mask level | [ ] |
| D71 | `mem -R -H -o -5 -v -T -h -z -a -C -Y -I` | usage_mem "Input/output" | `-a` outputs all secondaries; `-v` verbosity | [ ] |
| D72 | `mem -r -y -A -B -O -E -L -U -d -x -p -j -q -K -V -M -u` | usage_mem `*` options | accepted and ignored | [ ] |
| D73 | `mem` writes MD:Z | main_mem | `MB_F_WRITE_MD` always set | [ ] |

### Area D: names in minibwa.h

| ID | Name | Notes | Done |
|:---|:-----|:------|:-----|
| D80 | `MB_VERSION` | `--version`, `version`, @PG VN | [ ] |
| D81 | `MB_F_PAF` | `-f` | [ ] |
| D82 | `MB_F_NO_UNMAP` | `-u` | [ ] |
| D83 | `MB_F_COPY_COMMENT` | `-y` (A tests the option) | [ ] |
| D84 | `MB_F_PE` | default on with sr/adap (C) | [ ] |
| D85 | `MB_F_LONG` | B | [ ] |
| D86 | `MB_F_EQX` | B | [ ] |
| D87 | `MB_F_NO_KALLOC` | set by default unless built with `mimalloc=0` | [ ] |
| D88 | `MB_F_NO_ALN` | B | [ ] |
| D89 | `MB_F_PE_PREDEF` | C | [ ] |
| D90 | `MB_F_WRITE_DS`, `MB_F_WRITE_CS`, `MB_F_WRITE_MD` | B | [ ] |
| D91 | `MB_F_2ND_SEQ` | API only: SEQ/QUAL on secondary SAM records | [ ] |
| D92 | `MB_F_SUPP_SOFT` | `-Y` | [ ] |
| D93 | `MB_F_ADAP` | B | [ ] |
| D94 | `MB_F_PRIMARY5` | C | [ ] |
| D95 | `MB_F_NO_PAIRING` | C | [ ] |
| D96 | `MB_F_METH` | `--meth`; `mb_map_batch` reads it | [ ] |
| D97 | `MB_CIGAR_MATCH`, `MB_CIGAR_INS`, `MB_CIGAR_DEL`, `MB_CIGAR_N_SKIP`, `MB_CIGAR_SOFTCLIP`, `MB_CIGAR_HARDCLIP`, `MB_CIGAR_PADDING`, `MB_CIGAR_EQ_MATCH`, `MB_CIGAR_X_MISMATCH`, `MB_CIGAR_STR` | CIGAR decoding in the API; `=`/`X` with `--eqx` | [ ] |
| D98 | `mb_opt_t` | struct; field rows are A60-A66, B70-B86, C60-C64, D99 | [ ] |
| D99 | `mb_opt_t.flag`, `n_thread`, `out_n`, `out_s`, `seed`, `xa_max`, `cap_kalloc` | `seed` (11) feeds the tie-breaking hash; `n_thread` unused by the API | [ ] |
| D100 | `mb_idx_t` (opaque) | | [ ] |
| D101 | `mb_extra_t`: `cap`, `dp_score`, `dp_max0`, `dp_max`, `dp_max2`, `n_ambi`, `cs`, `n_cigar`, `cigar[]` | AS = dp_score, ms = dp_max0, md = dp_max - dp_max2; cs/MD string after `cigar[n_cigar]` when `cs` is set; `cap` internal | [ ] |
| D102 | `MB_PARENT_UNSET`, `MB_PARENT_TMP_PRI` | internal parent markers; not visible in results | n/a |
| D103 | `mb_hit_t`: `tid`, `ts`, `te`, `qs`, `qe`, `rev` | coordinates (query on the read's forward strand) | [ ] |
| D104 | `mb_hit_t`: `mapq`, `mlen`, `blen`, `p` | PAF columns 10-12 | [ ] |
| D105 | `mb_hit_t`: `parent`, `id`, `sam_pri` | primary iff parent == id; SAM primary iff sam_pri | [ ] |
| D106 | `mb_hit_t`: `score`, `score0`, `cnt`, `subsc`, `n_sub` | chaining score (s1), anchors (cm), s2 | [ ] |
| D107 | `mb_hit_t`: `inv`, `split`, `split_inv`, `frac_high`, `seed_ratio`, `rescued`, `proper_pair` | B (inv/split), C (rescued/proper_pair) | [ ] |
| D108 | `mb_hit_t`: `as`, `hash`, `flt`, `dummy` | "for internal uses"; do not pin | n/a |
| D109 | `mb_tbuf_t` (opaque) | | [ ] |
| D110 | `mb_idx_load`, `mb_idx_load_mmap`, `mb_idx_destroy`, `mb_idx_ctg_name`, `mb_idx_ctg_len` | rows D23-D26 | [ ] |
| D111 | `mb_opt_init`, `mb_opt_preset` | `mb_opt_init` = adap preset; B tests preset values | [ ] |
| D112 | `mb_tbuf_init`, `mb_tbuf_destroy`, `mb_tbuf_reset` | row D27 | [ ] |
| D113 | `mb_map`, `mb_map_batch` | rows D28, D29, C17 | [ ] |
