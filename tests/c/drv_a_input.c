/* Area A driver: reading queries (bseq.c), read-name matching (bseq.h),
 * the seeding and batching defaults of mb_opt_t, and the independence of
 * mb_map_batch() from how it cuts sub-batches (sb_seq, sb_len).
 *
 * The read cases print every record and batch as tab-separated lines for
 * pytest to compare with the input it wrote and with the batching rule. */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "mbpriv.h"
#include "minibwa_test.h"
/* bseq.o is part of the minibwa binary but not of libminibwa.a, which the
 * harness links drivers against, so compile the tree's bseq.c in here */
#include "bseq.c"

static void print_rec(const mb_bseq1_t *r)
{ // the comment goes last because it may contain tabs
	printf("R\t%llu\t%s\t%s\t%s\t%s\n", (unsigned long long)r->l_seq, r->name, r->seq, r->qual? r->qual : "*", r->comment? r->comment : "*");
}

static int check_rec(const mb_bseq1_t *r)
{
	int ok = 1;
	ok &= MBT_CHECKF(r->name && r->seq, "record without name or sequence");
	if (r->seq) ok &= MBT_CHECKF(strlen(r->seq) == r->l_seq, "%s: l_seq %llu, strlen %ld", r->name, (unsigned long long)r->l_seq, (long)strlen(r->seq));
	if (r->qual) ok &= MBT_CHECKF(strlen(r->qual) == r->l_seq, "%s: quality length differs", r->name);
	return ok;
}

static void free_rec(mb_bseq1_t *r)
{
	free(r->name); free(r->seq); free(r->qual); free(r->comment);
}

/* args: <file or -> <chunk_size> <with_qual> <with_comment> <frag_mode> <min_cnt> <max_chunk_size>
 * prints B batch n_records, then the records of the batch (R l_seq name seq qual comment) */
static int case_read(int argc, char *argv[])
{
	mb_bseq_file_t *fp;
	int64_t chunk, max_chunk;
	int with_qual, with_comment, frag, min_cnt, n, i, n_batch = 0;
	if (argc < 7) return MBT_FAIL("usage: read <file> <chunk> <with_qual> <with_comment> <frag> <min_cnt> <max_chunk>"), 1;
	chunk = atoll(argv[1]), with_qual = atoi(argv[2]), with_comment = atoi(argv[3]);
	frag = atoi(argv[4]), min_cnt = atoi(argv[5]), max_chunk = atoll(argv[6]);
	fp = mb_bseq_open(argv[0]);
	MBT_REQUIRE(fp != 0);
	while (1) {
		mb_bseq1_t *a;
		a = mb_bseq_read(fp, chunk, with_qual, with_comment, frag, min_cnt, max_chunk, &n);
		if (a == 0 || n == 0) {
			MBT_CHECK_INT(n, ==, 0);
			free(a);
			break;
		}
		printf("B\t%d\t%d\n", n_batch++, n);
		for (i = 0; i < n; ++i) {
			check_rec(&a[i]);
			print_rec(&a[i]);
			free_rec(&a[i]);
		}
		free(a);
	}
	mb_bseq_close(fp);
	MBT_CHECK(1);
	return 0;
}

/* args: <chunk_size> <with_qual> <with_comment> <file>...
 * reads the files in lockstep with mb_bseq_read_frag; prints B and R lines as `read` does */
static int case_frag(int argc, char *argv[])
{
	mb_bseq_file_t **fp;
	int64_t chunk;
	int with_qual, with_comment, n_fp, n, i, n_batch = 0;
	if (argc < 4) return MBT_FAIL("usage: frag <chunk> <with_qual> <with_comment> <file>..."), 1;
	chunk = atoll(argv[0]), with_qual = atoi(argv[1]), with_comment = atoi(argv[2]);
	n_fp = argc - 3;
	fp = (mb_bseq_file_t**)calloc(n_fp, sizeof(mb_bseq_file_t*));
	for (i = 0; i < n_fp; ++i) {
		fp[i] = mb_bseq_open(argv[3 + i]);
		MBT_REQUIRE(fp[i] != 0);
	}
	while (1) {
		mb_bseq1_t *a;
		a = mb_bseq_read_frag(n_fp, fp, chunk, with_qual, with_comment, &n);
		if (a == 0 || n == 0) {
			MBT_CHECK_INT(n, ==, 0);
			free(a);
			break;
		}
		MBT_CHECKF(n % n_fp == 0, "batch of %d records from %d files", n, n_fp);
		printf("B\t%d\t%d\n", n_batch++, n);
		for (i = 0; i < n; ++i) {
			check_rec(&a[i]);
			print_rec(&a[i]);
			free_rec(&a[i]);
		}
		free(a);
	}
	MBT_CHECK(mb_bseq_read_frag(0, fp, chunk, with_qual, with_comment, &n) == 0 && n == 0);
	for (i = 0; i < n_fp; ++i) mb_bseq_close(fp[i]);
	free(fp);
	return 0;
}

/* args: <path>...; prints O path opened(0/1) for mb_bseq_open */
static int case_open(int argc, char *argv[])
{
	int i;
	for (i = 0; i < argc; ++i) {
		mb_bseq_file_t *fp = mb_bseq_open(argv[i]);
		printf("O\t%s\t%d\n", argv[i], fp != 0);
		if (fp) mb_bseq_close(fp);
	}
	MBT_CHECK(1);
	return 0;
}

/* args: <name>...; prints Q name mb_qname_len, then S a b mb_qname_same for every ordered pair */
static int case_qname(int argc, char *argv[])
{
	int i, j;
	for (i = 0; i < argc; ++i) printf("Q\t%s\t%d\n", argv[i], mb_qname_len(argv[i]));
	for (i = 0; i < argc; ++i)
		for (j = 0; j < argc; ++j)
			printf("S\t%s\t%s\t%d\n", argv[i], argv[j], mb_qname_same(argv[i], argv[j]));
	MBT_CHECK(1);
	return 0;
}

/* the documented seeding and batching defaults (usage text: -k [19], -c [250], -K [100m,1g]) */
static int case_opt(int argc, char *argv[])
{
	mb_opt_t o;
	mb_opt_init(&o);
	MBT_CHECK_INT(o.min_len, ==, 19);
	MBT_CHECK_INT(o.max_occ, ==, 250);
	MBT_CHECK_INT(o.mb_size, ==, 100000000);
	MBT_CHECK_INT(o.max_mb_size, ==, 1000000000);
	MBT_CHECK(o.flag & MB_F_PE);
	MBT_CHECK(o.sb_seq > 0 && o.sb_len > 0);
	MBT_CHECK_INT(mb_opt_preset(&o, "sr"), ==, 0);
	MBT_CHECK_INT(o.mb_size, ==, 100000000);
	MBT_CHECK_INT(o.max_mb_size, ==, 1000000000);
	MBT_CHECK(o.flag & MB_F_PE);
	MBT_CHECK_INT(mb_opt_preset(&o, "lr"), ==, 0);
	MBT_CHECK_INT(o.mb_size, ==, 500000000);
	MBT_CHECK_INT(o.max_mb_size, ==, 1000000000);
	MBT_CHECK_INT(o.min_len, ==, 19);
	MBT_CHECK_INT(o.max_occ, ==, 250);
	MBT_CHECK(!(o.flag & MB_F_PE));
	return 0;
}

static int same_hits(int n1, const mb_hit_t *a, int n2, const mb_hit_t *b)
{
	int i, j;
	if (n1 != n2) return 0;
	for (i = 0; i < n1; ++i) {
		const mb_hit_t *x = &a[i], *y = &b[i];
		if (x->tid != y->tid || x->ts != y->ts || x->te != y->te || x->qs != y->qs || x->qe != y->qe || x->rev != y->rev)
			return 0;
		if (x->mapq != y->mapq || x->score != y->score || x->parent != y->parent || x->sam_pri != y->sam_pri)
			return 0;
		if ((x->p == 0) != (y->p == 0)) return 0;
		if (x->p) {
			if (x->p->n_cigar != y->p->n_cigar || x->p->dp_score != y->p->dp_score) return 0;
			for (j = 0; j < x->p->n_cigar; ++j)
				if (x->p->cigar[j] != y->p->cigar[j]) return 0;
		}
	}
	return 1;
}

static void free_hits(int n_seq, int32_t *n_hit, mb_hit_t **h)
{
	int i, j;
	for (i = 0; i < n_seq; ++i) {
		for (j = 0; j < n_hit[i]; ++j) free(h[i][j].p);
		free(h[i]);
	}
	free(h);
}

/* args: <idx_prefix> <reads.fa>
 * maps the reads single-end with mb_map_batch under several sub-batch limits (sb_seq,
 * sb_len) and checks that every read gets the same hits as with the defaults and as with
 * mb_map one read at a time */
static int case_subbatch(int argc, char *argv[])
{
	static const int32_t sb_seq[] = { 1, 2, 5, 24, 1000 };
	static const int32_t sb_len[] = { 1, 150, 1000, 1000000 };
	mb_opt_t opt;
	mb_idx_t *idx;
	mbt_seqs_t s;
	mb_hit_t **ref;
	int32_t *n_ref, *n_hit, i, k, l, n_mapped = 0;
	if (argc < 2) return MBT_FAIL("usage: subbatch <idx> <reads.fa>"), 1;
	idx = mb_idx_load(argv[0], 0);
	MBT_REQUIRE(idx != 0);
	MBT_REQUIRE(mbt_read_seqs(argv[1], &s) > 0);
	mb_opt_init(&opt);
	opt.flag &= ~MB_F_PE;
	n_ref = (int32_t*)calloc(s.n, sizeof(int32_t));
	n_hit = (int32_t*)calloc(s.n, sizeof(int32_t));
	ref = mb_map_batch(&opt, idx, s.n, s.len, (const char**)s.seq, n_ref, 0, (const char**)s.name);
	MBT_REQUIRE(ref != 0);
	for (i = 0; i < s.n; ++i) {
		int32_t n1;
		mb_hit_t *h1 = mb_map(&opt, idx, s.len[i], s.seq[i], 0, &n1, 0, s.name[i]);
		MBT_CHECKF(same_hits(n_ref[i], ref[i], n1, h1), "%s: mb_map and mb_map_batch differ", s.name[i]);
		for (k = 0; k < n1; ++k) free(h1[k].p);
		free(h1);
		n_mapped += n_ref[i] > 0;
	}
	MBT_CHECK(n_mapped > s.n / 2);
	for (k = 0; k < (int)(sizeof(sb_seq) / sizeof(sb_seq[0])); ++k)
		for (l = 0; l < (int)(sizeof(sb_len) / sizeof(sb_len[0])); ++l) {
			mb_opt_t o = opt;
			mb_hit_t **h;
			o.sb_seq = sb_seq[k], o.sb_len = sb_len[l];
			h = mb_map_batch(&o, idx, s.n, s.len, (const char**)s.seq, n_hit, 0, (const char**)s.name);
			MBT_REQUIRE(h != 0);
			for (i = 0; i < s.n; ++i)
				MBT_CHECKF(same_hits(n_ref[i], ref[i], n_hit[i], h[i]), "sb_seq=%d sb_len=%d: %s differs", o.sb_seq, o.sb_len, s.name[i]);
			free_hits(s.n, n_hit, h);
		}
	free_hits(s.n, n_ref, ref);
	free(n_hit);
	mbt_free_seqs(&s);
	mb_idx_destroy(idx);
	return 0;
}

static const mbt_case_t cases[] = {
	{ "read", case_read },
	{ "frag", case_frag },
	{ "open", case_open },
	{ "qname", case_qname },
	{ "opt", case_opt },
	{ "subbatch", case_subbatch },
	{ 0, 0 }
};

MBT_MAIN(cases)
