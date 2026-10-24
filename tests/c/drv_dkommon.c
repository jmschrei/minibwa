/* Area D driver: the helpers every output byte goes through (kommon.c,
 * mb_escape), the bisulfite helpers of l2bit.c, and kthread.c (not part of
 * libminibwa.a; the harness links kthread.o from libmbt_extra.a).
 *
 * Cases:
 *   sprintf_lite   kom_sprintf_lite/km_sprintf_lite against snprintf
 *   parse_num      kom_parse_num with k/m/g suffixes and decimals
 *   escape         mb_escape on \t \n \r \\ and a trailing backslash
 *   strdup         kom_strdup
 *   meth_convert   l2b_meth_convert, l2b_meth_rev
 *   intv2cid_meth  l2b_intv2cid_meth on the four converted copies of a --meth index
 *   kt_for         every index processed once, thread ids in range
 *   kt_pipeline    steps run in order and outputs come out in input order
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <limits.h>
#include <pthread.h>
#include "kommon.h"
#include "kalloc.h"
#include "l2bit.h"
#include "mbpriv.h"
#include "kthread.h"
#include "minibwa_test.h"

static int case_sprintf_lite(int argc, char *argv[])
{
	static const int ints[] = { 0, 1, -1, 9, 10, -10, 99, 100, 12345, -12345, 2147483647, -2147483647, 1000000000 };
	static const long longs[] = { 0L, 1L, -1L, 4294967296L, -4294967296L, 9223372036854775807L, -9223372036854775807L, 1234567890123L };
	static const uint32_t uints[] = { 0U, 1U, 15U, 16U, 4294967295U, 2147483648U, 3000000000U };
	static const char *strs[] = { "", "a", "chr1", "a much longer string that forces the buffer to grow more than once ................................................" };
	kstring_t s = {0, 0, 0}, t = {0, 0, 0};
	char want[4096], *big;
	size_t i, j, total = 0;
	int64_t len;
	void *km;
	(void)argc, (void)argv;
	for (i = 0; i < sizeof(ints) / sizeof(ints[0]); ++i) {
		s.l = 0;
		len = kom_sprintf_lite(&s, "x%dy", ints[i]);
		snprintf(want, sizeof(want), "x%dy", ints[i]);
		MBT_CHECK_STR(s.s, want);
		MBT_CHECK_INT(len, ==, (long long)strlen(want));
		MBT_CHECK_INT(s.l, ==, strlen(want));
	}
	for (i = 0; i < sizeof(longs) / sizeof(longs[0]); ++i) {
		s.l = 0;
		kom_sprintf_lite(&s, "%ld", longs[i]);
		snprintf(want, sizeof(want), "%ld", longs[i]);
		MBT_CHECK_STR(s.s, want);
	}
	for (i = 0; i < sizeof(uints) / sizeof(uints[0]); ++i) {
		s.l = 0;
		kom_sprintf_lite(&s, "%u,", uints[i]);
		snprintf(want, sizeof(want), "%u,", (unsigned)uints[i]);
		MBT_CHECK_STR(s.s, want);
	}
	for (i = 0; i < sizeof(strs) / sizeof(strs[0]); ++i) {
		s.l = 0;
		len = kom_sprintf_lite(&s, "[%s]%c%s", strs[i], 'Z', strs[i]);
		snprintf(want, sizeof(want), "[%s]%c%s", strs[i], 'Z', strs[i]);
		MBT_CHECK_STR(s.s, want);
		MBT_CHECK_INT(len, ==, (long long)strlen(want));
	}
	// a SAM-like line; appending keeps what is there; NUL termination after every call
	s.l = 0;
	kom_sprintf_lite(&s, "%s\t%d\t%s\t%ld\t%d\t", "read1", 2064, "chr2", 123456789L, 60);
	kom_sprintf_lite(&s, "%dM%c", 150, 'H');
	kom_sprintf_lite(&s, "\tNM:i:%d\tXA:Z:%s,%c%ld,%u;", 3, "chrX", '-', 77L, 4000000000U);
	snprintf(want, sizeof(want), "%s\t%d\t%s\t%ld\t%d\t%dM%c\tNM:i:%d\tXA:Z:%s,%c%ld,%u;", "read1", 2064, "chr2", 123456789L, 60, 150, 'H', 3,
		"chrX", '-', 77L, 4000000000U);
	MBT_CHECK_STR(s.s, want);
	MBT_CHECK_INT(s.l, ==, strlen(want));
	// literal text only (the NULL-buffer mode that only returns a length has no caller, so it is not tested)
	s.l = 0;
	MBT_CHECK_INT(kom_sprintf_lite(&s, "no conversions"), ==, 14);
	MBT_CHECK_STR(s.s, "no conversions");
	// many appends grow the buffer; compare with one snprintf of the whole
	s.l = 0;
	big = (char*)malloc(1 << 20);
	for (i = 0, j = 0; i < 20000; ++i) {
		kom_sprintf_lite(&s, "%d,%s;", (int)(i * 7919) - 50000, (i & 1)? "ab" : "");
		j += snprintf(big + j, (1 << 20) - j, "%d,%s;", (int)(i * 7919) - 50000, (i & 1)? "ab" : "");
		total = j;
	}
	MBT_CHECK_INT(s.l, ==, total);
	MBT_CHECK(s.l == total && memcmp(s.s, big, total) == 0 && s.s[s.l] == 0);
	// km_sprintf_lite with and without a kalloc pool
	km = km_init();
	for (i = 0; i < 2; ++i) {
		void *k = i == 0? 0 : km;
		t.s = 0, t.l = t.m = 0;
		km_sprintf_lite(k, &t, "%s:%d-%ld", "ctg", 5, 6000000000L);
		MBT_CHECK_STR(t.s, "ctg:5-6000000000");
		if (k == 0) free(t.s); // with a pool, the string is in the pool (HAVE_KALLOC) or leaked here
	}
	km_destroy(km);
	free(big);
	free(s.s);
	return 0;
}

static int case_parse_num(int argc, char *argv[])
{
	static const struct { const char *s; int64_t v; const char *rest; } t[] = {
		{ "0", 0, "" }, { "10", 10, "" }, { "1k", 1000, "" }, { "1K", 1000, "" }, { "1.5k", 1500, "" }, { "2m", 2000000, "" },
		{ "2M", 2000000, "" }, { "1g", 1000000000, "" }, { "1G", 1000000000, "" }, { "0.25G", 250000000, "" }, { "100m,1g", 100000000, ",1g" },
		{ "0.6", 1, "" }, { "0.4", 0, "" }, { "2.6", 3, "" }, { "2.4", 2, "" }, { "7x", 7, "x" }, { "0.001k", 1, "" }, { "5k,", 5000, "," }
	};
	size_t i;
	(void)argc, (void)argv;
	for (i = 0; i < sizeof(t) / sizeof(t[0]); ++i) {
		char *rest = 0;
		int64_t v = kom_parse_num(t[i].s, &rest);
		MBT_CHECKF(v == t[i].v, "kom_parse_num(\"%s\") = %lld, want %lld", t[i].s, (long long)v, (long long)t[i].v);
		MBT_CHECKF(rest && strcmp(rest, t[i].rest) == 0, "kom_parse_num(\"%s\") stops at \"%s\", want \"%s\"", t[i].s, rest? rest : "(null)", t[i].rest);
		MBT_CHECK_INT(kom_parse_num(t[i].s, 0), ==, t[i].v);
	}
	return 0;
}

static int case_escape(int argc, char *argv[])
{
	static const struct { const char *in, *out; } t[] = {
		{ "", "" }, { "plain", "plain" }, { "@RG\\tID:x", "@RG\tID:x" }, { "a\\nb\\rc", "a\nb\rc" }, { "back\\\\slash", "back\\slash" },
		{ "\\t\\t", "\t\t" }, { "end\\", "end" }, { "\\\\t", "\\t" }
	};
	size_t i;
	(void)argc, (void)argv;
	for (i = 0; i < sizeof(t) / sizeof(t[0]); ++i) {
		char *s = kom_strdup(t[i].in);
		MBT_CHECK(mb_escape(s) == s);
		MBT_CHECK_STR(s, t[i].out);
		free(s);
	}
	return 0;
}

static int case_strdup(int argc, char *argv[])
{
	const char *src = "a string\twith a tab";
	char *d = kom_strdup(src);
	(void)argc, (void)argv;
	MBT_CHECK(d != src);
	MBT_CHECK_STR(d, src);
	free(d);
	d = kom_strdup("");
	MBT_CHECK_STR(d, "");
	free(d);
	return 0;
}

static int case_meth_convert(int argc, char *argv[])
{
	uint8_t seq[10] = { 0, 1, 2, 3, 4, 1, 1, 2, 2, 0 }, s2[10];
	static const uint8_t c2t[10] = { 0, 3, 2, 3, 4, 3, 3, 2, 2, 0 }, g2a[10] = { 0, 1, 0, 3, 4, 1, 1, 0, 0, 0 };
	(void)argc, (void)argv;
	memcpy(s2, seq, 10); l2b_meth_convert(L2B_METH_C2T, 10, s2);
	MBT_CHECK(memcmp(s2, c2t, 10) == 0);
	memcpy(s2, seq, 10); l2b_meth_convert(L2B_METH_G2A, 10, s2);
	MBT_CHECK(memcmp(s2, g2a, 10) == 0);
	memcpy(s2, seq, 10); l2b_meth_convert(L2B_METH_NONE, 10, s2);
	MBT_CHECK(memcmp(s2, seq, 10) == 0);
	memcpy(s2, seq, 10); l2b_meth_convert(L2B_METH_C2T, 3, s2); // only the first len bases
	MBT_CHECK(s2[1] == 3 && s2[5] == 1);
	MBT_CHECK(l2b_meth_rev(L2B_METH_NONE) == L2B_METH_NONE);
	MBT_CHECK(l2b_meth_rev(L2B_METH_C2T) == L2B_METH_G2A);
	MBT_CHECK(l2b_meth_rev(L2B_METH_G2A) == L2B_METH_C2T);
	return 0;
}

/* args: <meth_index_prefix> <n_ctg> <len1> ...
 * The --meth BWT concatenates, each tot_len long: C->T forward, G->A forward, the reverse
 * complement of G->A, the reverse complement of C->T (l2b_save_pac_meth). */
static int case_intv2cid_meth(int argc, char *argv[])
{
	char fn[4096];
	l2b_t *l2b;
	int32_t n_ctg, i, k;
	int64_t T, off;
	if (argc < 2) return MBT_FAIL("usage: intv2cid_meth <prefix> <n_ctg> <len>..."), 1;
	snprintf(fn, sizeof(fn), "%s.l2b", argv[0]);
	MBT_REQUIRE((l2b = l2b_load(fn)) != 0);
	n_ctg = atoi(argv[1]);
	MBT_REQUIRE(argc == 2 + n_ctg && (int32_t)l2b->n_ctg == n_ctg);
	for (i = 0, T = 0; i < n_ctg; ++i) T += atol(argv[2 + i]);
	MBT_REQUIRE((int64_t)l2b->tot_len == T);
	for (i = 0, off = 0; i < n_ctg; off += atol(argv[2 + i]), ++i) {
		int64_t L = atol(argv[2 + i]);
		int64_t pos_list[3] = { 0, L / 2, L - 1 };
		for (k = 0; k < 3; ++k) {
			int64_t pos = pos_list[k], len = pos + 25 <= L? 25 : L - pos;
			int64_t st[4], cst = -1;
			int32_t c, rev = -1;
			l2b_meth_t mt;
			st[0] = off + pos;
			st[1] = T + off + pos;
			st[2] = 3 * T - (off + pos + len);
			st[3] = 4 * T - (off + pos + len);
			for (c = 0; c < 4; ++c) {
				int64_t tid = l2b_intv2cid_meth(l2b, st[c], st[c] + len, &mt, &cst, &rev);
				MBT_CHECKF(tid == i && cst == pos && rev == (c >= 2) && mt == (c == 0 || c == 3? L2B_METH_C2T : L2B_METH_G2A),
					"contig %d pos %ld len %ld copy %d: tid %ld cst %ld rev %d mt %d", i, (long)pos, (long)len, c, (long)tid, (long)cst, rev, (int)mt);
			}
		}
		if (i + 1 < n_ctg) { // an interval across the end of this contig is in no contig
			int64_t cst, e = off + L;
			int rev;
			l2b_meth_t mt;
			MBT_CHECK(l2b_intv2cid_meth(l2b, e - 5, e + 5, &mt, &cst, &rev) < 0);
			MBT_CHECK(l2b_intv2cid_meth(l2b, 3 * T - e - 5, 3 * T - e + 5, &mt, &cst, &rev) < 0);
		}
	}
	{
		int64_t cst;
		int rev;
		l2b_meth_t mt;
		MBT_CHECK(l2b_intv2cid_meth(l2b, 4 * T - 3, 4 * T + 2, &mt, &cst, &rev) < 0); // past the end
	}
	l2b_destroy(l2b);
	return 0;
}

typedef struct {
	long n;
	int n_threads, bad_tid;
	int *count;
	pthread_mutex_t lock;
} for_data_t;

static void for_func(void *data, long i, int tid)
{
	for_data_t *d = (for_data_t*)data;
	pthread_mutex_lock(&d->lock);
	if (i >= 0 && i < d->n) ++d->count[i];
	if (tid < 0 || tid >= d->n_threads) d->bad_tid = 1;
	pthread_mutex_unlock(&d->lock);
}

static int case_kt_for(int argc, char *argv[])
{
	static const int threads[] = { 1, 2, 3, 8 };
	static const long ns[] = { 0, 1, 2, 7, 1000, 10007 };
	size_t a, b;
	(void)argc, (void)argv;
	for (a = 0; a < sizeof(threads) / sizeof(threads[0]); ++a) {
		for (b = 0; b < sizeof(ns) / sizeof(ns[0]); ++b) {
			for_data_t d;
			long i, bad = 0;
			memset(&d, 0, sizeof(d));
			d.n = ns[b], d.n_threads = threads[a];
			d.count = (int*)calloc(d.n + 1, sizeof(int));
			pthread_mutex_init(&d.lock, 0);
			kt_for(d.n_threads, for_func, &d, d.n);
			for (i = 0; i < d.n; ++i)
				if (d.count[i] != 1) ++bad;
			MBT_CHECKF(bad == 0 && !d.bad_tid, "kt_for(%d threads, n=%ld): %ld indices not processed exactly once", d.n_threads, d.n, bad);
			pthread_mutex_destroy(&d.lock);
			free(d.count);
		}
	}
	return 0;
}

typedef struct {
	int n_items, next_in, n_out, err;
	int *out;
} pl_shared_t;

typedef struct {
	int id, steps_done;
} pl_item_t;

static void *pl_func(void *shared, int step, void *in)
{
	pl_shared_t *p = (pl_shared_t*)shared;
	if (step == 0) {
		pl_item_t *it;
		if (p->next_in >= p->n_items) return 0;
		it = (pl_item_t*)calloc(1, sizeof(pl_item_t));
		it->id = p->next_in++;
		it->steps_done = 1;
		return it;
	} else if (step == 1) {
		pl_item_t *it = (pl_item_t*)in;
		volatile long x = 0, k;
		for (k = 0; k < (it->id % 5) * 2000; ++k) x += k; // uneven work
		if (it->steps_done != 1) p->err = 1;
		it->steps_done = 2;
		return it;
	} else {
		pl_item_t *it = (pl_item_t*)in;
		if (it->steps_done != 2) p->err = 1;
		p->out[p->n_out++] = it->id;
		free(it);
		return 0;
	}
}

static int case_kt_pipeline(int argc, char *argv[])
{
	static const int threads[] = { 1, 2, 3, 4 };
	size_t a;
	(void)argc, (void)argv;
	for (a = 0; a < sizeof(threads) / sizeof(threads[0]); ++a) {
		pl_shared_t p;
		int i, bad = 0;
		memset(&p, 0, sizeof(p));
		p.n_items = 300;
		p.out = (int*)calloc(p.n_items, sizeof(int));
		kt_pipeline(threads[a], pl_func, &p, 3);
		MBT_CHECK_INT(p.n_out, ==, p.n_items);
		for (i = 0; i < p.n_out; ++i)
			if (p.out[i] != i) ++bad;
		MBT_CHECKF(bad == 0 && !p.err, "kt_pipeline(%d threads): %d items out of order, step order error %d", threads[a], bad, p.err);
		free(p.out);
	}
	return 0;
}

static const mbt_case_t cases[] = {
	{ "sprintf_lite", case_sprintf_lite },
	{ "parse_num", case_parse_num },
	{ "escape", case_escape },
	{ "strdup", case_strdup },
	{ "meth_convert", case_meth_convert },
	{ "intv2cid_meth", case_intv2cid_meth },
	{ "kt_for", case_kt_for },
	{ "kt_pipeline", case_kt_pipeline },
	{ 0, 0 }
};

MBT_MAIN(cases)
