/* Area A driver: BWT, suffix-array, l2bit and seeding lookups.
 *
 * Each case prints what it computed as tab-separated lines for pytest to
 * compare with a naive suffix array of the indexed text (tests/helpers_a.py),
 * and checks here only what is defined by another call of minibwa: rank
 * functions against each other, batched against unbatched SA and SMEM
 * lookups, mb_seed_intv_batch against mb_seed_intv, and mb_anchor_sort's
 * order. The k-mer cache is checked only through the SMEMs found with it:
 * its table layout is not part of any contract.
 * Intervals are compared as sets: the order in which batched calls return
 * them is not part of their contract. */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "mbpriv.h"
#include "minibwa_test.h"

static uint64_t drv_rng = 11;

static uint64_t drv_rand(void)
{
	return kom_splitmix64(&drv_rng);
}

static uint8_t **encode_seqs(const mbt_seqs_t *s)
{
	int32_t i, j;
	uint8_t **q;
	q = (uint8_t**)calloc(s->n, sizeof(uint8_t*));
	for (i = 0; i < s->n; ++i) {
		q[i] = (uint8_t*)malloc(s->len[i] + 1);
		for (j = 0; j < s->len[i]; ++j)
			q[i][j] = kom_nt4_table[(uint8_t)s->seq[i][j]];
	}
	return q;
}

static void free_seqs4(int32_t n, uint8_t **q)
{
	int32_t i;
	for (i = 0; i < n; ++i) free(q[i]);
	free(q);
}

static int cmp_sai(const void *a_, const void *b_)
{
	const mb_sai_t *a = (const mb_sai_t*)a_, *b = (const mb_sai_t*)b_;
	if (a->info != b->info) return a->info < b->info? -1 : 1;
	if (a->x[0] != b->x[0]) return a->x[0] < b->x[0]? -1 : 1;
	if (a->x[1] != b->x[1]) return a->x[1] < b->x[1]? -1 : 1;
	if (a->size != b->size) return a->size < b->size? -1 : 1;
	return 0;
}

static size_t sort_uniq_sai(size_t n, mb_sai_t *a)
{
	size_t i, j;
	if (n == 0) return 0;
	qsort(a, n, sizeof(mb_sai_t), cmp_sai);
	for (i = 1, j = 0; i < n; ++i)
		if (cmp_sai(&a[i], &a[j]) != 0) a[++j] = a[i];
	return j + 1;
}

static int same_sai_set(size_t n1, mb_sai_t *a1, size_t n2, mb_sai_t *a2)
{
	size_t i;
	n1 = sort_uniq_sai(n1, a1);
	n2 = sort_uniq_sai(n2, a2);
	if (n1 != n2) return 0;
	for (i = 0; i < n1; ++i)
		if (cmp_sai(&a1[i], &a2[i]) != 0) return 0;
	return 1;
}

static void print_sai(const char *tag, const char *name, const mb_sai_t *p)
{
	printf("%s\t%s\t%d\t%d\t%llu\t%llu\t%llu\n", tag, name, (int32_t)(p->info>>32), (int32_t)p->info,
		(unsigned long long)p->x[0], (unsigned long long)p->x[1], (unsigned long long)p->size);
}

/* args: <in.mbw>
 * prints N seq_len primary L2[1..4], then R k cnt[0..3] = mb_bwt_rank1a(k) for every row
 * boundary k in [0, seq_len+1]; checks mb_bwt_rank11 and mb_bwt_rank2a against it. The pair
 * (primary, primary+1) is skipped: base minibwa's rank2a counts the base of the row after
 * the sentinel there, which makes forward extension run past the end of the text */
static int case_rank(int argc, char *argv[])
{
	mb_bwt_t *bwt;
	uint64_t k, l, n, (*cnt)[4], tk[4], tl[4];
	int c, i;
	if (argc < 1) return MBT_FAIL("usage: rank <in.mbw>"), 1;
	bwt = mb_bwt_load(argv[0]);
	MBT_REQUIRE(bwt != 0);
	n = bwt->seq_len;
	printf("N\t%llu\t%llu", (unsigned long long)n, (unsigned long long)bwt->primary);
	for (c = 1; c <= 4; ++c) printf("\t%llu", (unsigned long long)bwt->L2[c]);
	putchar('\n');
	cnt = (uint64_t(*)[4])calloc(n + 2, sizeof(*cnt));
	for (k = 0; k <= n + 1; ++k) {
		mb_bwt_rank1a(bwt, k, cnt[k]);
		printf("R\t%llu\t%llu\t%llu\t%llu\t%llu\n", (unsigned long long)k, (unsigned long long)cnt[k][0],
			(unsigned long long)cnt[k][1], (unsigned long long)cnt[k][2], (unsigned long long)cnt[k][3]);
	}
	for (k = 0; k <= n + 1; ++k)
		for (c = 0; c < 4; ++c) {
			uint64_t r = mb_bwt_rank11(bwt, k, c);
			if (r != cnt[k][c]) MBT_FAIL("rank11(%llu, %d) = %llu, rank1a gives %llu", (unsigned long long)k, c, (unsigned long long)r, (unsigned long long)cnt[k][c]);
		}
	for (k = 0; k <= n + 1; ++k)
		for (l = k; l <= n + 1 && l <= k + 300; ++l) {
			if (k == bwt->primary && l == k + 1) continue; // wrong in base minibwa: counts the base after the sentinel row
			mb_bwt_rank2a(bwt, k, l, tk, tl);
			for (c = 0; c < 4; ++c)
				if (tk[c] != cnt[k][c] || tl[c] != cnt[l][c])
					break;
			if (c < 4) MBT_FAIL("rank2a(%llu, %llu) differs from rank1a", (unsigned long long)k, (unsigned long long)l);
		}
	for (i = 0; i < 20000; ++i) {
		k = drv_rand() % (n + 2), l = drv_rand() % (n + 2);
		if (k > l) { uint64_t t = k; k = l; l = t; }
		if (k == bwt->primary && l == k + 1) continue;
		mb_bwt_rank2a(bwt, k, l, tk, tl);
		for (c = 0; c < 4; ++c)
			if (tk[c] != cnt[k][c] || tl[c] != cnt[l][c])
				break;
		if (c < 4) MBT_FAIL("rank2a(%llu, %llu) differs from rank1a", (unsigned long long)k, (unsigned long long)l);
	}
	MBT_CHECK(1);
	free(cnt);
	mb_bwt_destroy(bwt);
	return 0;
}

/* args: <in.mbw> <batch_size>...
 * prints S row sa = mb_bwt_sa(row) for every row; checks mb_bwt_sa_batch in consecutive
 * chunks of each batch size and on a shuffled list with repeats */
static int case_sa(int argc, char *argv[])
{
	mb_bwt_t *bwt;
	uint64_t k, n, *s, *x;
	int i, j, n_bad = 0;
	if (argc < 1) return MBT_FAIL("usage: sa <in.mbw> [batch_size...]"), 1;
	bwt = mb_bwt_load(argv[0]);
	MBT_REQUIRE(bwt != 0);
	n = bwt->seq_len;
	s = (uint64_t*)calloc(n + 1, sizeof(uint64_t));
	x = (uint64_t*)calloc(n + 1, sizeof(uint64_t));
	for (k = 0; k <= n; ++k) {
		s[k] = mb_bwt_sa(bwt, k);
		printf("S\t%llu\t%lld\n", (unsigned long long)k, (long long)s[k]);
	}
	for (i = 1; i < argc; ++i) {
		int64_t b = atol(argv[i]), st, m;
		MBT_REQUIRE(b > 0);
		for (st = 0; st <= (int64_t)n; st += b) {
			m = st + b <= (int64_t)n + 1? b : (int64_t)n + 1 - st;
			for (j = 0; j < m; ++j) x[j] = st + j;
			mb_bwt_sa_batch(0, bwt, m, x);
			for (j = 0; j < m; ++j)
				if (x[j] != s[st + j] && n_bad++ < 10)
					MBT_FAIL("batch %lld: sa_batch(row %lld) = %lld, sa() = %lld", (long long)b, (long long)(st + j), (long long)x[j], (long long)s[st + j]);
		}
	}
	for (i = 0; i < 50; ++i) { // shuffled rows with repeats
		uint64_t rows[61];
		int m = 1 + drv_rand() % 61;
		for (j = 0; j < m; ++j) x[j] = rows[j] = j > 0 && drv_rand() % 5 == 0? rows[j - 1] : drv_rand() % (n + 1);
		mb_bwt_sa_batch(0, bwt, m, x);
		for (j = 0; j < m; ++j)
			if (x[j] != s[rows[j]] && n_bad++ < 10)
				MBT_FAIL("shuffled sa_batch(row %llu) = %lld, sa() = %lld", (unsigned long long)rows[j], (long long)x[j], (long long)s[rows[j]]);
	}
	mb_bwt_sa_batch(0, bwt, 0, x); // empty batch is a no-op
	MBT_CHECK_INT(n_bad, ==, 0);
	free(x); free(s);
	mb_bwt_destroy(bwt);
	return 0;
}

static void backward_intv(const mb_bwt_t *bwt, int len, const uint8_t *p, mb_sai_t *ik)
{
	mb_sai_t ok[4];
	int i;
	mb_bwt_set_intv(bwt, p[len - 1], ik);
	for (i = len - 2; i >= 0; --i) {
		mb_bwt_extend(bwt, ik, ok, 1);
		*ik = ok[p[i]];
	}
}

static void forward_intv(const mb_bwt_t *bwt, int len, const uint8_t *p, mb_sai_t *ik)
{
	mb_sai_t ok[4];
	int i;
	mb_bwt_set_intv(bwt, p[0], ik);
	for (i = 1; i < len; ++i) {
		mb_bwt_extend(bwt, ik, ok, 0);
		*ik = ok[3 - p[i]];
	}
}

/* args: <in.mbw> <max_len>
 * prints X pattern x0 x1 size for every pattern of 1..max_len bases, built by backward
 * extension; checks that forward extension gives the same interval */
static int case_extend(int argc, char *argv[])
{
	mb_bwt_t *bwt;
	int len, max_len, i;
	uint8_t p[16];
	char str[17];
	if (argc < 2) return MBT_FAIL("usage: extend <in.mbw> <max_len>"), 1;
	bwt = mb_bwt_load(argv[0]);
	MBT_REQUIRE(bwt != 0);
	max_len = atoi(argv[1]);
	MBT_REQUIRE(max_len >= 1 && max_len <= 8);
	for (len = 1; len <= max_len; ++len) {
		uint64_t z, nz = 1ULL << 2 * len;
		for (z = 0; z < nz; ++z) {
			mb_sai_t b, f;
			for (i = 0; i < len; ++i) p[i] = z >> 2 * (len - 1 - i) & 3, str[i] = "ACGT"[p[i]];
			str[len] = 0;
			backward_intv(bwt, len, p, &b);
			forward_intv(bwt, len, p, &f);
			if (b.size != f.size || (b.size > 0 && (b.x[0] != f.x[0] || b.x[1] != f.x[1])))
				MBT_FAIL("%s: backward [%llu,%llu,%llu] forward [%llu,%llu,%llu]", str, (unsigned long long)b.x[0],
					(unsigned long long)b.x[1], (unsigned long long)b.size, (unsigned long long)f.x[0], (unsigned long long)f.x[1], (unsigned long long)f.size);
			printf("X\t%s\t%llu\t%llu\t%llu\n", str, (unsigned long long)b.x[0], (unsigned long long)b.x[1], (unsigned long long)b.size);
		}
	}
	MBT_CHECK(1);
	mb_bwt_destroy(bwt);
	return 0;
}

static size_t smem_one(const mb_bwt_t *bwt, int32_t st, int32_t en, const uint8_t *q, int32_t min_len, int32_t min_occ, mb_sai_t **a, size_t *m)
{
	int64_t x = st;
	size_t n = 0;
	mb_sai_t p;
	if (st >= en) return 0;
	do {
		x = mb_bwt_smem(bwt, en, q, x, min_len, min_occ, &p);
		if (p.size > 0) {
			kom_grow(mb_sai_t, *a, n, *m);
			(*a)[n++] = p;
		}
	} while (x < en);
	return n;
}

/* args: <in.mbw> <reads.fa> <min_len> <min_occ> <cache_len> <batch_size>...
 * prints M name qs qe x0 x1 size for the SMEMs mb_bwt_smem finds in each read, and
 * m name ... for those in the middle half [len/4, len-len/4); checks that
 * mb_bwt_smem_batch finds the same sets for each batch size, with whole reads and with
 * whole reads and middle halves interleaved in one batch */
static int case_smem(int argc, char *argv[])
{
	mb_bwt_t *bwt;
	mbt_seqs_t s;
	uint8_t **q;
	int32_t i, j, k, min_len, min_occ, cache_len;
	mb_sai_t **one, **mid;
	size_t *n_one, *n_mid, m;
	if (argc < 5) return MBT_FAIL("usage: smem <in.mbw> <reads.fa> <min_len> <min_occ> <cache_len> [batch_size...]"), 1;
	bwt = mb_bwt_load(argv[0]);
	MBT_REQUIRE(bwt != 0);
	MBT_REQUIRE(mbt_read_seqs(argv[1], &s) > 0);
	min_len = atoi(argv[2]), min_occ = atoi(argv[3]), cache_len = atoi(argv[4]);
	if (cache_len > 0) mb_bwt_cache(bwt, cache_len);
	q = encode_seqs(&s);
	one = (mb_sai_t**)calloc(s.n, sizeof(mb_sai_t*));
	mid = (mb_sai_t**)calloc(s.n, sizeof(mb_sai_t*));
	n_one = (size_t*)calloc(s.n, sizeof(size_t));
	n_mid = (size_t*)calloc(s.n, sizeof(size_t));
	for (i = 0; i < s.n; ++i) {
		int32_t st = s.len[i] / 4, en = s.len[i] - s.len[i] / 4;
		m = 0;
		n_one[i] = smem_one(bwt, 0, s.len[i], q[i], min_len, min_occ, &one[i], &m);
		for (j = 0; j < (int32_t)n_one[i]; ++j) print_sai("M", s.name[i], &one[i][j]);
		m = 0;
		n_mid[i] = smem_one(bwt, st, en, q[i], min_len, min_occ, &mid[i], &m);
		for (j = 0; j < (int32_t)n_mid[i]; ++j) print_sai("m", s.name[i], &mid[i][j]);
	}
	for (k = 5; k < argc; ++k) {
		int32_t b = atoi(argv[k]), mode;
		mb_smem_entry_t *e;
		mb_sai_v *v;
		MBT_REQUIRE(b > 0);
		e = (mb_smem_entry_t*)calloc(b, sizeof(mb_smem_entry_t));
		v = (mb_sai_v*)calloc(b, sizeof(mb_sai_v));
		for (mode = 0; mode < 2; ++mode) { // 0: whole reads; 1: whole reads and middle halves alternating
			int32_t n_ent = mode == 0? s.n : 2 * s.n;
			for (i = 0; i < n_ent; i += b) {
				int32_t n = i + b <= n_ent? b : n_ent - i;
				memset(e, 0, b * sizeof(mb_smem_entry_t));
				for (j = 0; j < n; ++j) {
					int32_t r = mode == 0? i + j : (i + j) / 2, half = mode == 1 && (i + j) % 2 == 1;
					e[j].min_len = min_len, e[j].min_occ = min_occ;
					e[j].st = half? s.len[r] / 4 : 0;
					e[j].en = half? s.len[r] - s.len[r] / 4 : s.len[r];
					e[j].q = q[r];
					e[j].v = &v[j];
					v[j].n = 0;
				}
				mb_bwt_smem_batch(0, bwt, n, e);
				for (j = 0; j < n; ++j) {
					int32_t r = mode == 0? i + j : (i + j) / 2, half = mode == 1 && (i + j) % 2 == 1;
					mb_sai_t *want = half? mid[r] : one[r];
					size_t n_want = half? n_mid[r] : n_one[r];
					if (half && e[j].st >= e[j].en) n_want = 0;
					if (!same_sai_set(n_want, want, v[j].n, v[j].a))
						MBT_FAIL("batch size %d: mb_bwt_smem_batch finds %ld intervals in %s%s, mb_bwt_smem %ld",
							b, (long)v[j].n, s.name[r], half? " (middle half)" : "", (long)n_want);
				}
			}
		}
		for (j = 0; j < b; ++j) free(v[j].a);
		free(v); free(e);
	}
	MBT_CHECK(1);
	for (i = 0; i < s.n; ++i) free(one[i]), free(mid[i]);
	free(one); free(mid); free(n_one); free(n_mid);
	free_seqs4(s.n, q);
	mbt_free_seqs(&s);
	mb_bwt_destroy(bwt);
	return 0;
}

/* args: <idx_prefix> <reads.fa> <min_len> <max_sub_occ> <group_size>...
 * loads the index as minibwa map does (with the 10-mer cache); prints D name qs qe x0 x1 size
 * for the set of intervals mb_seed_intv returns for each read; checks that
 * mb_seed_intv_batch returns the same set for reads processed in groups of each size */
static int case_seed(int argc, char *argv[])
{
	mb_idx_t *idx;
	mbt_seqs_t s;
	uint8_t **q;
	int32_t i, j, k, min_len, max_sub_occ;
	mb_sai_v *one;
	if (argc < 4) return MBT_FAIL("usage: seed <idx> <reads.fa> <min_len> <max_sub_occ> [group_size...]"), 1;
	idx = mb_idx_load(argv[0], 0);
	MBT_REQUIRE(idx != 0);
	MBT_REQUIRE(mbt_read_seqs(argv[1], &s) > 0);
	min_len = atoi(argv[2]), max_sub_occ = atoi(argv[3]);
	q = encode_seqs(&s);
	one = (mb_sai_v*)calloc(s.n, sizeof(mb_sai_v));
	for (i = 0; i < s.n; ++i) {
		mb_seed_intv(0, idx->bwt, s.len[i], q[i], min_len, max_sub_occ, &one[i]);
		one[i].n = sort_uniq_sai(one[i].n, one[i].a);
		for (j = 0; j < (int32_t)one[i].n; ++j) print_sai("D", s.name[i], &one[i].a[j]);
	}
	for (k = 4; k < argc; ++k) {
		int32_t g = atoi(argv[k]);
		MBT_REQUIRE(g > 0);
		for (i = 0; i < s.n; i += g) {
			int32_t n = i + g <= s.n? g : s.n - i;
			mb_sai_v *v = (mb_sai_v*)calloc(n, sizeof(mb_sai_v));
			mb_seed_intv_batch(0, idx->bwt, n, &s.len[i], &q[i], min_len, max_sub_occ, v);
			for (j = 0; j < n; ++j) {
				if (!same_sai_set(one[i + j].n, one[i + j].a, v[j].n, v[j].a))
					MBT_FAIL("group size %d: mb_seed_intv_batch gives %ld intervals for %s, mb_seed_intv %ld", g, (long)v[j].n, s.name[i + j], (long)one[i + j].n);
				free(v[j].a);
			}
			free(v);
		}
	}
	MBT_CHECK(1);
	for (i = 0; i < s.n; ++i) free(one[i].a);
	free(one);
	free_seqs4(s.n, q);
	mbt_free_seqs(&s);
	mb_idx_destroy(idx);
	return 0;
}

static int cmp_anchor_key(const mb_anchor_t *a, const mb_anchor_t *b)
{ // the order mb_anchor and mb_anchor_sort promise: contig, then strand, then target end
	if (a->sid != b->sid) return a->sid < b->sid? -1 : 1;
	if (a->tpos != b->tpos) return a->tpos < b->tpos? -1 : 1;
	return 0;
}

static int cmp_anchor_full(const void *a_, const void *b_)
{
	const mb_anchor_t *a = (const mb_anchor_t*)a_, *b = (const mb_anchor_t*)b_;
	int c = cmp_anchor_key(a, b);
	if (c != 0) return c;
	if (a->qpos != b->qpos) return a->qpos < b->qpos? -1 : 1;
	if (a->len != b->len) return a->len < b->len? -1 : 1;
	return 0;
}

/* args: <idx_prefix> <reads.fa> <min_len> <max_sub_occ> <max_occ>
 * prints A name tid strand qs ts len for the anchors mb_anchor derives from mb_seed_intv's
 * intervals (qs on the read as given, ts on the forward strand, both 0-based); checks the
 * return value (1.0 without --meth), the order, and that mb_anchor_sort restores that
 * order from a shuffled copy */
static int case_anchor(int argc, char *argv[])
{
	mb_idx_t *idx;
	mbt_seqs_t s;
	uint8_t **q;
	int32_t i, min_len, max_sub_occ, max_occ;
	if (argc < 5) return MBT_FAIL("usage: anchor <idx> <reads.fa> <min_len> <max_sub_occ> <max_occ>"), 1;
	idx = mb_idx_load(argv[0], 0);
	MBT_REQUIRE(idx != 0);
	MBT_REQUIRE(mbt_read_seqs(argv[1], &s) > 0);
	min_len = atoi(argv[2]), max_sub_occ = atoi(argv[3]), max_occ = atoi(argv[4]);
	q = encode_seqs(&s);
	for (i = 0; i < s.n; ++i) {
		mb_sai_v u = {0,0,0};
		mb_anchor_v v = {0,0,0};
		mb_anchor_t *b, *c;
		double ratio;
		int64_t j;
		int32_t qlen = s.len[i];
		mb_seed_intv(0, idx->bwt, qlen, q[i], min_len, max_sub_occ, &u);
		ratio = mb_anchor(0, idx, &u, min_len, qlen, q[i], L2B_METH_NONE, max_occ, &v);
		MBT_CHECKF(ratio == 1.0, "%s: mb_anchor returned %g", s.name[i], ratio);
		for (j = 0; j < v.n; ++j) {
			const mb_anchor_t *a = &v.a[j];
			int32_t rev = a->sid & 1;
			int32_t qs = rev? qlen - 1 - a->qpos : a->qpos + 1 - a->len;
			printf("A\t%s\t%d\t%c\t%d\t%lld\t%d\n", s.name[i], a->sid >> 1, "+-"[rev], qs, (long long)(a->tpos + 1 - a->len), a->len);
			if (j > 0 && cmp_anchor_key(&v.a[j - 1], a) > 0)
				MBT_FAIL("%s: anchors %ld and %ld are out of order", s.name[i], (long)j - 1, (long)j);
		}
		if (v.n > 0) { // mb_anchor_sort on a shuffled copy
			b = (mb_anchor_t*)malloc(v.n * sizeof(mb_anchor_t));
			c = (mb_anchor_t*)malloc(v.n * sizeof(mb_anchor_t));
			memcpy(b, v.a, v.n * sizeof(mb_anchor_t));
			for (j = v.n - 1; j > 0; --j) {
				int64_t r = drv_rand() % (j + 1);
				mb_anchor_t t = b[j]; b[j] = b[r]; b[r] = t;
			}
			mb_anchor_sort(idx->l2b, v.n, b);
			for (j = 1; j < v.n; ++j)
				if (cmp_anchor_key(&b[j - 1], &b[j]) > 0) break;
			MBT_CHECKF(j >= v.n, "%s: mb_anchor_sort left anchor %ld out of order", s.name[i], (long)j);
			memcpy(c, v.a, v.n * sizeof(mb_anchor_t));
			qsort(b, v.n, sizeof(mb_anchor_t), cmp_anchor_full);
			qsort(c, v.n, sizeof(mb_anchor_t), cmp_anchor_full);
			for (j = 0; j < v.n; ++j)
				if (cmp_anchor_full(&b[j], &c[j]) != 0) break;
			MBT_CHECKF(j >= v.n, "%s: mb_anchor_sort changed the anchors", s.name[i]);
			free(b); free(c);
		}
		free(u.a); free(v.a);
	}
	MBT_CHECK(1);
	free_seqs4(s.n, q);
	mbt_free_seqs(&s);
	mb_idx_destroy(idx);
	return 0;
}

/* args: <in.l2b> <len>...
 * prints the contig table (H, C), the ambiguous intervals (B), every contig through
 * l2b_getseq (G), l2b_getambi on a grid of windows (W) and l2b_intv2cid for every start
 * on both strands and each length (I) */
static int case_l2b(int argc, char *argv[])
{
	l2b_t *l2b;
	uint64_t i;
	int64_t tid, st, en, ret, tot2;
	int a, n_bad = 0;
	if (argc < 1) return MBT_FAIL("usage: l2b <in.l2b> [len...]"), 1;
	l2b = l2b_load(argv[0]);
	MBT_REQUIRE(l2b != 0);
	printf("H\t%llu\t%llu\t%llu\n", (unsigned long long)l2b->n_ctg, (unsigned long long)l2b->tot_len, (unsigned long long)l2b->n_ambi);
	for (i = 0; i < l2b->n_ctg; ++i)
		printf("C\t%llu\t%s\t%llu\t%llu\t%s\n", (unsigned long long)i, l2b->ctg[i].name, (unsigned long long)l2b->ctg[i].len,
			(unsigned long long)l2b->ctg[i].off, l2b->ctg[i].comm? l2b->ctg[i].comm : "*");
	for (i = 0; i < l2b->n_ambi; ++i)
		printf("B\t%llu\t%llu\n", (unsigned long long)l2b->ambi[i].st, (unsigned long long)l2b->ambi[i].en);
	for (tid = 0; tid < (int64_t)l2b->n_ctg; ++tid) {
		int64_t len = l2b->ctg[tid].len, j;
		uint8_t *full = (uint8_t*)malloc(len + 1), *buf = (uint8_t*)malloc(len + 1);
		ret = l2b_getseq(l2b, tid, 0, len, full);
		MBT_CHECK_INT(ret, ==, len);
		printf("G\t%lld\t", (long long)tid);
		for (j = 0; j < len; ++j) putchar('0' + full[j]);
		putchar('\n');
		ret = l2b_getseq(l2b, tid, -5, len + 5, buf); // clamped to the contig
		MBT_CHECK_INT(ret, ==, len);
		MBT_CHECK(len == 0 || memcmp(buf, full, len) == 0);
		for (st = 0; st <= len; st += len < 80? 1 : 13) // every sub-window agrees with the whole contig
			for (en = st; en <= len + 3; en += len < 80? 1 : 17) {
				int64_t e = en < len? en : len;
				ret = l2b_getseq(l2b, tid, st, en, buf);
				if ((ret != e - st || (e > st && memcmp(buf, full + st, e - st) != 0)) && n_bad++ < 10)
					MBT_FAIL("l2b_getseq(%lld, %lld, %lld) returned %lld or differs from the contig", (long long)tid, (long long)st, (long long)en, (long long)ret);
				if (en < len) {
					int32_t n_ambi;
					ret = l2b_getambi(l2b, tid, st, en, &n_ambi);
					printf("W\t%lld\t%lld\t%lld\t%lld\t%d\n", (long long)tid, (long long)st, (long long)en, (long long)ret, n_ambi);
				}
			}
		free(full); free(buf);
	}
	MBT_CHECK_INT(n_bad, ==, 0);
	MBT_CHECK_INT(l2b_getseq(l2b, -1, 0, 1, 0), ==, -1);
	MBT_CHECK_INT(l2b_getseq(l2b, l2b->n_ctg, 0, 1, 0), ==, -1);
	{
		int32_t n_ambi = 7;
		MBT_CHECK_INT(l2b_getambi(l2b, -1, 0, 1, &n_ambi), ==, -1);
		MBT_CHECK_INT(n_ambi, ==, 0);
		n_ambi = 7;
		MBT_CHECK_INT(l2b_getambi(l2b, l2b->n_ctg, 0, 1, &n_ambi), ==, -1);
		MBT_CHECK_INT(n_ambi, ==, 0);
	}
	tot2 = l2b->tot_len * 2;
	for (a = 1; a < argc; ++a) {
		int64_t len = atol(argv[a]);
		for (st = 0; st <= tot2; ++st) {
			int64_t cst = -1;
			int rev = -1;
			ret = l2b_intv2cid(l2b, st, st + len, &cst, &rev);
			if (ret >= 0) printf("I\t%lld\t%lld\t%lld\t%lld\t%d\n", (long long)st, (long long)(st + len), (long long)ret, (long long)cst, rev);
			else printf("I\t%lld\t%lld\t%lld\t.\t.\n", (long long)st, (long long)(st + len), (long long)ret);
		}
	}
	l2b_destroy(l2b);
	return 0;
}

/* prints T byte code for every byte of kom_nt4_table */
static int case_nt4(int argc, char *argv[])
{
	int i;
	for (i = 0; i < 256; ++i) printf("T\t%d\t%d\n", i, kom_nt4_table[i]);
	MBT_CHECK(1);
	return 0;
}

/* args: <in.mbw>; mb_bwt_load returns NULL for a missing file or a wrong magic number */
static int case_load(int argc, char *argv[])
{
	mb_bwt_t *bwt;
	l2b_t *l2b;
	if (argc < 3) return MBT_FAIL("usage: load <good.mbw> <good.l2b> <bad_file>"), 1;
	bwt = mb_bwt_load(argv[0]);
	MBT_REQUIRE(bwt != 0);
	MBT_CHECK_INT(bwt->L2[4], ==, bwt->seq_len);
	mb_bwt_destroy(bwt);
	mb_bwt_destroy(0);
	l2b = l2b_load(argv[1]);
	MBT_REQUIRE(l2b != 0);
	l2b_destroy(l2b);
	MBT_CHECK(mb_bwt_load(argv[2]) == 0);
	MBT_CHECK(l2b_load(argv[2]) == 0);
	MBT_CHECK(mb_bwt_load(argv[1]) == 0); // an .l2b is not an .mbw
	MBT_CHECK(l2b_load(argv[0]) == 0);
	MBT_CHECK(mb_bwt_load("/nonexistent/x.mbw") == 0);
	MBT_CHECK(l2b_load("/nonexistent/x.l2b") == 0);
	MBT_CHECK(mb_idx_load("/nonexistent/x", 0) == 0);
	return 0;
}

static const mbt_case_t cases[] = {
	{ "rank", case_rank },
	{ "sa", case_sa },
	{ "extend", case_extend },
	{ "smem", case_smem },
	{ "seed", case_seed },
	{ "anchor", case_anchor },
	{ "l2b", case_l2b },
	{ "nt4", case_nt4 },
	{ "load", case_load },
	{ 0, 0 }
};

MBT_MAIN(cases)
