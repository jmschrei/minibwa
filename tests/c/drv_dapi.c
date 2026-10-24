/* Area D driver: the libminibwa API (minibwa.h) and the SAM/PAF writer it
 * shares with `minibwa map`.
 *
 * Cases:
 *   idx_load     mb_idx_load/mb_idx_load_mmap/mb_idx_destroy/mb_idx_ctg_*
 *   opt          mb_opt_init/mb_opt_preset
 *   dump         map reads through mb_map() or mb_map_batch() in one of several
 *                modes and print every hit in a PAF-like layout for pytest to
 *                compare with `minibwa map -f`
 *   fmt          map reads with mb_map() and print the SAM/PAF records that
 *                mb_format() writes for them under given output flags
 *   batch_edge   mb_map_batch() with n_seq <= 0 or a meth flag on a plain index
 *   seed         opt.seed feeds the hash that breaks ties between equal hits
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "minibwa.h"
#include "mbpriv.h"
#include "minibwa_test.h"

/* check the documented relations between the fields of one read's hits */
static void check_hits(const mb_idx_t *idx, const char *name, int32_t qlen, int32_t n_hit, const mb_hit_t *hit)
{
	int32_t j, n_sam_pri = 0;
	for (j = 0; j < n_hit; ++j) {
		const mb_hit_t *h = &hit[j];
		int64_t ctg_len = mb_idx_ctg_len(idx, h->tid);
		int32_t k, ql = 0, tl = 0;
		MBT_CHECKF(mb_idx_ctg_name(idx, h->tid) != 0 && ctg_len > 0, "%s: hit %d tid %ld out of range", name, j, (long)h->tid);
		MBT_CHECKF(0 <= h->ts && h->ts < h->te && h->te <= ctg_len, "%s: hit %d spans %ld-%ld", name, j, (long)h->ts, (long)h->te);
		MBT_CHECKF(0 <= h->qs && h->qs < h->qe && h->qe <= qlen, "%s: hit %d query %d-%d of %d", name, j, h->qs, h->qe, qlen);
		MBT_CHECKF(h->id == j, "%s: hit %d has id %d", name, j, h->id);
		MBT_CHECKF(h->parent >= 0 && h->parent < n_hit, "%s: hit %d parent %d", name, j, h->parent);
		if (h->parent >= 0 && h->parent < n_hit)
			MBT_CHECKF(hit[h->parent].parent == h->parent, "%s: hit %d's parent %d is not a primary", name, j, h->parent);
		MBT_CHECKF(0 <= h->mlen && h->mlen <= h->blen, "%s: hit %d mlen %d blen %d", name, j, h->mlen, h->blen);
		MBT_CHECKF(0 <= h->mapq && h->mapq <= 60, "%s: hit %d mapq %d", name, j, h->mapq);
		if (h->sam_pri) {
			++n_sam_pri;
			MBT_CHECKF(h->parent == h->id, "%s: SAM primary hit %d is a secondary", name, j);
		}
		if (!MBT_CHECKF(h->p != 0, "%s: hit %d has no base alignment", name, j)) continue;
		for (k = 0; k < h->p->n_cigar; ++k) {
			int op = h->p->cigar[k] & 0xf, len = h->p->cigar[k] >> 4;
			MBT_CHECKF(op == MB_CIGAR_MATCH || op == MB_CIGAR_INS || op == MB_CIGAR_DEL || op == MB_CIGAR_EQ_MATCH || op == MB_CIGAR_X_MISMATCH,
				"%s: hit %d CIGAR op %d", name, j, op);
			if (op == MB_CIGAR_MATCH || op == MB_CIGAR_INS || op == MB_CIGAR_EQ_MATCH || op == MB_CIGAR_X_MISMATCH) ql += len;
			if (op == MB_CIGAR_MATCH || op == MB_CIGAR_DEL || op == MB_CIGAR_EQ_MATCH || op == MB_CIGAR_X_MISMATCH) tl += len;
		}
		MBT_CHECKF(ql == h->qe - h->qs && tl == h->te - h->ts, "%s: hit %d CIGAR spans %d/%d, hit %d/%ld", name, j, ql, tl,
			h->qe - h->qs, (long)(h->te - h->ts));
	}
	MBT_CHECKF(n_hit == 0 || n_sam_pri == 1, "%s: %d SAM-primary hits among %d", name, n_sam_pri, n_hit);
}

/* one line per hit: the fields of `minibwa map -f` plus sam_pri */
static void print_hits(const mb_idx_t *idx, const char *name, int32_t qlen, int32_t n_hit, const mb_hit_t *hit)
{
	int32_t j, k;
	printf("N\t%s\t%d\n", name, n_hit);
	for (j = 0; j < n_hit; ++j) {
		const mb_hit_t *h = &hit[j];
		int pri = h->parent == h->id;
		printf("H\t%s\t%d\t%d\t%d\t%c\t%s\t%ld\t%ld\t%ld\t%d\t%d\t%d\ttp:A:%c\ts1:i:%d\tcm:i:%d", name, qlen, h->qs, h->qe, "+-"[h->rev],
			mb_idx_ctg_name(idx, h->tid), (long)mb_idx_ctg_len(idx, h->tid), (long)h->ts, (long)h->te, h->mlen, h->blen, h->mapq,
			pri? 'P' : 'S', h->score, h->cnt);
		if (pri) printf("\ts2:i:%d", h->subsc >= 0? h->subsc : 0);
		if (h->p) {
			printf("\tNM:i:%d\tAS:i:%d\tms:i:%d\tmd:i:%d\tcg:Z:", h->blen - h->mlen + h->p->n_ambi, h->p->dp_score, h->p->dp_max0,
				h->p->dp_max - h->p->dp_max2);
			for (k = 0; k < h->p->n_cigar; ++k)
				printf("%d%c", h->p->cigar[k] >> 4, MB_CIGAR_STR[h->p->cigar[k] & 0xf]);
			if (h->p->cs) printf("\t%s", (const char*)&h->p->cigar[h->p->n_cigar]);
		}
		printf("\tsp:i:%d\n", h->sam_pri);
	}
}

static void free_hits(int32_t n_hit, mb_hit_t *hit)
{
	int32_t j;
	for (j = 0; j < n_hit; ++j) free(hit[j].p);
	free(hit);
}

static uint64_t parse_flag(const char *s)
{
	if (strcmp(s, "eqx") == 0) return MB_F_EQX;
	if (strcmp(s, "cs") == 0) return MB_F_WRITE_CS;
	if (strcmp(s, "ds") == 0) return MB_F_WRITE_DS;
	if (strcmp(s, "MD") == 0) return MB_F_WRITE_MD;
	if (strcmp(s, "paf") == 0) return MB_F_PAF;
	if (strcmp(s, "2nd") == 0) return MB_F_2ND_SEQ;
	if (strcmp(s, "soft") == 0) return MB_F_SUPP_SOFT;
	if (strcmp(s, "comment") == 0) return MB_F_COPY_COMMENT;
	if (strcmp(s, "meth") == 0) return MB_F_METH;
	if (strcmp(s, "nopair") == 0) return MB_F_NO_PAIRING;
	if (strcmp(s, "pe") == 0) return MB_F_PE;
	return 0;
}

/* trailing arguments: flag names, and key=value for mt, seed, nthread, sbseq, sblen, mmap, outs, xa */
static int parse_opts(int argc, char *argv[], mb_opt_t *opt, int *mt, int *mmap)
{
	int i;
	for (i = 0; i < argc; ++i) {
		if (strncmp(argv[i], "mt=", 3) == 0) *mt = atoi(argv[i] + 3);
		else if (strncmp(argv[i], "mmap=", 5) == 0) *mmap = atoi(argv[i] + 5);
		else if (strncmp(argv[i], "outs=", 5) == 0) opt->out_s = atof(argv[i] + 5);
		else if (strncmp(argv[i], "xa=", 3) == 0) opt->xa_max = atoi(argv[i] + 3);
		else if (strncmp(argv[i], "seed=", 5) == 0) opt->seed = atoi(argv[i] + 5);
		else if (strncmp(argv[i], "nthread=", 8) == 0) opt->n_thread = atoi(argv[i] + 8);
		else if (strncmp(argv[i], "sbseq=", 6) == 0) opt->sb_seq = atoi(argv[i] + 6);
		else if (strncmp(argv[i], "sblen=", 6) == 0) opt->sb_len = atoi(argv[i] + 6);
		else if (strcmp(argv[i], "-pe") == 0) opt->flag &= ~MB_F_PE;
		else if (parse_flag(argv[i]) != 0) opt->flag |= parse_flag(argv[i]);
		else return MBT_FAIL("unknown option %s", argv[i]), -1;
	}
	return 0;
}

/* args: <prefix> <missing_prefix> <meth_prefix> <n_ctg> <name1> <len1> ... */
static int case_idx_load(int argc, char *argv[])
{
	mb_idx_t *idx, *idx2;
	int32_t i, n_ctg, k;
	if (argc < 4) return MBT_FAIL("usage: idx_load <prefix> <missing> <meth_prefix> <n_ctg> [name len]..."), 1;
	n_ctg = atoi(argv[3]);
	MBT_REQUIRE(argc == 4 + 2 * n_ctg);
	MBT_CHECK(mb_idx_load(argv[1], 0) == 0);
	MBT_CHECK(mb_idx_load(argv[1], 1) == 0);
	MBT_CHECK(mb_idx_load_mmap(argv[1], 0, 1) == 0);
	MBT_CHECK(mb_idx_load_mmap(argv[1], 0, 0) == 0);
	MBT_CHECK(mb_idx_load(argv[0], 1) == 0); // no .meth.mbw next to a plain index
	MBT_CHECK(mb_idx_load_mmap(argv[0], 1, 1) == 0);
	mb_idx_destroy(0); // NULL-safe
	for (k = 0; k < 5; ++k) {
		const char *pfx = k < 3? argv[0] : argv[2];
		idx = k == 0? mb_idx_load(pfx, 0) : k == 1? mb_idx_load_mmap(pfx, 0, 1) : k == 2? mb_idx_load_mmap(pfx, 0, 0)
			: k == 3? mb_idx_load(pfx, 1) : mb_idx_load_mmap(pfx, 1, 1);
		if (!MBT_CHECKF(idx != 0, "loader %d failed on %s", k, pfx)) continue;
		for (i = 0; i < n_ctg; ++i) {
			MBT_CHECK_STR(mb_idx_ctg_name(idx, i), argv[4 + 2 * i]);
			MBT_CHECK_INT(mb_idx_ctg_len(idx, i), ==, atol(argv[5 + 2 * i]));
		}
		MBT_CHECK(mb_idx_ctg_name(idx, -1) == 0);
		MBT_CHECK(mb_idx_ctg_name(idx, n_ctg) == 0);
		MBT_CHECK(mb_idx_ctg_name(idx, n_ctg + 1000) == 0);
		MBT_CHECK_INT(mb_idx_ctg_len(idx, -1), ==, -1);
		MBT_CHECK_INT(mb_idx_ctg_len(idx, n_ctg), ==, -1);
		mb_idx_destroy(idx);
	}
	// the meth index also loads as a plain index
	idx2 = mb_idx_load(argv[2], 0);
	MBT_CHECK(idx2 != 0);
	mb_idx_destroy(idx2);
	return 0;
}

static int case_opt(int argc, char *argv[])
{
	mb_opt_t a, b, c;
	(void)argc, (void)argv;
	mb_opt_init(&a);
	memset(&b, 0xff, sizeof(b)); // mb_opt_preset() must reset every field
	MBT_CHECK_INT(mb_opt_preset(&b, "adap"), ==, 0);
	MBT_CHECK(memcmp(&a, &b, sizeof(mb_opt_t)) == 0);
	MBT_CHECK((a.flag & (MB_F_PE|MB_F_ADAP)) == (MB_F_PE|MB_F_ADAP));
	MBT_CHECK((a.flag & (MB_F_PAF|MB_F_NO_UNMAP|MB_F_COPY_COMMENT|MB_F_LONG|MB_F_EQX|MB_F_NO_ALN|MB_F_PE_PREDEF|MB_F_WRITE_DS|MB_F_WRITE_CS|MB_F_WRITE_MD|MB_F_2ND_SEQ|MB_F_SUPP_SOFT|MB_F_PRIMARY5|MB_F_NO_PAIRING|MB_F_METH)) == 0);
	// defaults stated in the usage text of `minibwa map`
	MBT_CHECK_INT(a.n_thread, ==, 1);
	MBT_CHECK_INT(a.out_n, ==, 0);
	MBT_CHECK(a.out_s == 0.8f);
	MBT_CHECK_INT(a.xa_max, ==, 5);
	MBT_CHECK_INT(a.seed, ==, 11);
	MBT_CHECK_INT(a.mb_size, ==, 100000000);
	MBT_CHECK_INT(a.max_mb_size, ==, 1000000000);
	MBT_CHECK_INT(mb_opt_preset(&c, "sr"), ==, 0);
	MBT_CHECK((c.flag & MB_F_PE) && !(c.flag & MB_F_ADAP));
	MBT_CHECK_INT(mb_opt_preset(&c, "lr"), ==, 0);
	MBT_CHECK((c.flag & MB_F_LONG) && !(c.flag & MB_F_PE) && !(c.flag & MB_F_ADAP));
	MBT_CHECK_INT(mb_opt_preset(&c, "nonexistent"), ==, -1);
	MBT_CHECK_INT(mb_opt_preset(&c, ""), ==, -1);
	return 0;
}

/* args: <prefix> <is_meth> <mode> <reads.fa> [opts...]
 * modes: one, one_tbuf, one_tbuf_nok, one_reset, one_nt4, one_parity (mt 1/2 by read parity),
 *        batch, batch_tbuf */
static int case_dump(int argc, char *argv[])
{
	mb_opt_t opt;
	mb_idx_t *idx;
	mbt_seqs_t s;
	mb_tbuf_t *b = 0;
	const char *mode;
	int32_t i, mt = 0, use_mmap = -1;
	if (argc < 4) return MBT_FAIL("usage: dump <prefix> <is_meth> <mode> <reads.fa> [opts]"), 1;
	mode = argv[2];
	mb_opt_init(&opt);
	if (parse_opts(argc - 4, argv + 4, &opt, &mt, &use_mmap) < 0) return 1;
	idx = use_mmap < 0? mb_idx_load(argv[0], atoi(argv[1])) : mb_idx_load_mmap(argv[0], atoi(argv[1]), use_mmap);
	MBT_REQUIRE(idx != 0);
	MBT_REQUIRE(mbt_read_seqs(argv[3], &s) > 0);
	if (strcmp(mode, "one_tbuf") == 0 || strcmp(mode, "one_reset") == 0 || strcmp(mode, "batch_tbuf") == 0) b = mb_tbuf_init(0);
	else if (strcmp(mode, "one_tbuf_nok") == 0) b = mb_tbuf_init(1);
	if (strncmp(mode, "one", 3) == 0) {
		for (i = 0; i < s.n; ++i) {
			mb_hit_t *hit;
			int32_t n_hit = -1, mti = strcmp(mode, "one_parity") == 0? 1 + (i&1) : mt;
			char *seq = s.seq[i];
			if (strcmp(mode, "one_nt4") == 0) { // 0/1/2/3 encoding; anything else as 4
				int32_t l;
				seq = (char*)malloc(s.len[i]);
				for (l = 0; l < s.len[i]; ++l) {
					const char *p = strchr("ACGTacgt", s.seq[i][l]);
					seq[l] = p && s.seq[i][l]? (char)((p - "ACGTacgt") & 3) : 4;
				}
			}
			hit = mb_map(&opt, idx, s.len[i], seq, mti, &n_hit, b, s.name[i]);
			MBT_CHECKF(n_hit >= 0 && (n_hit == 0 || hit != 0), "%s: n_hit %d, hit %p", s.name[i], n_hit, (void*)hit);
			check_hits(idx, s.name[i], s.len[i], n_hit, hit);
			print_hits(idx, s.name[i], s.len[i], n_hit, hit);
			free_hits(n_hit, hit);
			if (seq != s.seq[i]) free(seq);
			if (b && strcmp(mode, "one_reset") == 0) mb_tbuf_reset(b, 0);
			else if (b) mb_tbuf_reset(b, opt.cap_kalloc);
		}
	} else if (strncmp(mode, "batch", 5) == 0) {
		mb_hit_t **hit;
		int32_t *n_hit = (int32_t*)calloc(s.n, sizeof(int32_t));
		hit = mb_map_batch(&opt, idx, s.n, s.len, (const char**)s.seq, n_hit, b, (const char**)s.name);
		MBT_REQUIRE(hit != 0);
		for (i = 0; i < s.n; ++i) {
			check_hits(idx, s.name[i], s.len[i], n_hit[i], hit[i]);
			print_hits(idx, s.name[i], s.len[i], n_hit[i], hit[i]);
			free_hits(n_hit[i], hit[i]);
		}
		free(hit);
		free(n_hit);
	} else return MBT_FAIL("unknown mode %s", mode), 1;
	if (b) mb_tbuf_destroy(b);
	mbt_free_seqs(&s);
	mb_idx_destroy(idx);
	return 0;
}

/* args: <prefix> <reads.fq> [opts...]: SAM (or PAF with "paf") from mb_format() for every hit,
 * after the filtering `minibwa map` applies with --outn 1000000 --outs 0 (none). */
static int case_fmt(int argc, char *argv[])
{
	mb_opt_t opt;
	mb_idx_t *idx;
	mbt_seqs_t s;
	FILE *fp;
	int32_t i, mt = 0, use_mmap = -1;
	kstring_t out = {0, 0, 0};
	char *line = 0;
	size_t cap = 0;
	if (argc < 2) return MBT_FAIL("usage: fmt <prefix> <reads.fq> [opts]"), 1;
	mb_opt_init(&opt);
	if (parse_opts(argc - 2, argv + 2, &opt, &mt, &use_mmap) < 0) return 1;
	idx = mb_idx_load(argv[0], 0);
	MBT_REQUIRE(idx != 0);
	MBT_REQUIRE(mbt_read_seqs(argv[1], &s) > 0);
	// quality strings: the fourth line of each FASTQ record
	MBT_REQUIRE((fp = fopen(argv[1], "r")) != 0);
	for (i = 0; i < s.n; ++i) {
		mb_bseq1_t t;
		mb_hit_t *hit, *hv[1];
		int32_t n_hit, j, k;
		ssize_t l;
		for (k = 0; k < 4; ++k)
			MBT_REQUIRE((l = getline(&line, &cap, fp)) > 0);
		if (line[l - 1] == '\n') line[--l] = 0;
		MBT_REQUIRE(l == s.len[i]);
		memset(&t, 0, sizeof(t));
		t.name = s.name[i], t.seq = s.seq[i], t.qual = line, t.l_seq = s.len[i];
		hit = mb_map(&opt, idx, s.len[i], s.seq[i], 0, &n_hit, 0, s.name[i]);
		hv[0] = hit;
		out.l = 0;
		if (n_hit == 0) mb_format(0, &out, idx->l2b, &t, 1, &n_hit, hv, -1, &opt, 0, 0);
		for (j = 0; j < n_hit; ++j)
			mb_format(0, &out, idx->l2b, &t, 1, &n_hit, hv, j, &opt, 0, 0);
		fwrite(out.s, 1, out.l, stdout);
		free_hits(n_hit, hit);
	}
	free(line);
	free(out.s);
	fclose(fp);
	mbt_free_seqs(&s);
	mb_idx_destroy(idx);
	return 0;
}

/* args: <plain_prefix> <reads.fa> */
static int case_batch_edge(int argc, char *argv[])
{
	mb_opt_t opt;
	mb_idx_t *idx;
	mbt_seqs_t s;
	int32_t n_hit[4] = {7, 7, 7, 7}, n1 = -1;
	mb_hit_t **hit, *h1;
	if (argc < 2) return MBT_FAIL("usage: batch_edge <prefix> <reads.fa>"), 1;
	mb_opt_init(&opt);
	opt.flag &= ~MB_F_PE;
	idx = mb_idx_load(argv[0], 0);
	MBT_REQUIRE(idx != 0);
	MBT_REQUIRE(mbt_read_seqs(argv[1], &s) > 0);
	MBT_CHECK(mb_map_batch(&opt, idx, 0, s.len, (const char**)s.seq, n_hit, 0, (const char**)s.name) == 0);
	MBT_CHECK(mb_map_batch(&opt, idx, -1, s.len, (const char**)s.seq, n_hit, 0, (const char**)s.name) == 0);
	opt.flag |= MB_F_METH; // a meth flag needs an index loaded with is_meth
	MBT_CHECK(mb_map_batch(&opt, idx, 1, s.len, (const char**)s.seq, n_hit, 0, (const char**)s.name) == 0);
	opt.flag &= ~MB_F_METH;
	// mt != 0 on a plain index: no hits
	h1 = mb_map(&opt, idx, s.len[0], s.seq[0], 1, &n1, 0, s.name[0]);
	MBT_CHECK(h1 == 0 && n1 == 0);
	n1 = -1;
	h1 = mb_map(&opt, idx, s.len[0], s.seq[0], 2, &n1, 0, s.name[0]);
	MBT_CHECK(h1 == 0 && n1 == 0);
	// qname may be NULL in mb_map_batch
	hit = mb_map_batch(&opt, idx, 1, s.len, (const char**)s.seq, n_hit, 0, 0);
	MBT_REQUIRE(hit != 0);
	MBT_CHECK(n_hit[0] > 0);
	check_hits(idx, s.name[0], s.len[0], n_hit[0], hit[0]);
	free_hits(n_hit[0], hit[0]);
	free(hit);
	mbt_free_seqs(&s);
	mb_idx_destroy(idx);
	return 0;
}

/* args: <prefix> <reads.fa> <n_seeds>: print, per seed value, the SAM-primary position of each read */
static int case_seed(int argc, char *argv[])
{
	mb_opt_t opt;
	mb_idx_t *idx;
	mbt_seqs_t s;
	int32_t i, j, sd, n_seeds;
	if (argc < 3) return MBT_FAIL("usage: seed <prefix> <reads.fa> <n_seeds>"), 1;
	n_seeds = atoi(argv[2]);
	idx = mb_idx_load(argv[0], 0);
	MBT_REQUIRE(idx != 0);
	MBT_REQUIRE(mbt_read_seqs(argv[1], &s) > 0);
	for (sd = 0; sd < n_seeds; ++sd) {
		mb_opt_init(&opt);
		opt.seed = 11 + sd * 1000003;
		printf("S\t%d", opt.seed);
		for (i = 0; i < s.n; ++i) {
			int32_t n_hit;
			mb_hit_t *hit = mb_map(&opt, idx, s.len[i], s.seq[i], 0, &n_hit, 0, s.name[i]);
			check_hits(idx, s.name[i], s.len[i], n_hit, hit);
			for (j = 0; j < n_hit; ++j)
				if (hit[j].sam_pri) printf("\t%s:%ld:%c", mb_idx_ctg_name(idx, hit[j].tid), (long)hit[j].ts, "+-"[hit[j].rev]);
			free_hits(n_hit, hit);
		}
		putchar('\n');
	}
	mbt_free_seqs(&s);
	mb_idx_destroy(idx);
	return 0;
}

static const mbt_case_t cases[] = {
	{ "idx_load", case_idx_load },
	{ "opt", case_opt },
	{ "dump", case_dump },
	{ "fmt", case_fmt },
	{ "batch_edge", case_batch_edge },
	{ "seed", case_seed },
	{ 0, 0 }
};

MBT_MAIN(cases)
