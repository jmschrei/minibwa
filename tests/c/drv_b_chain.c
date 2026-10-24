/* Area B driver: chaining (mb_lchain_dp) and the chain-to-hit conversions in
 * map-algo.c (mb_cal_mblen, mb_gen_hit, mb_split_hit, mb_squeeze_a,
 * mb_sync_hits), on hand-built anchors with hand-worked expected values, and
 * mb_lchain_dp against a brute-force chaining DP on random anchors.
 *
 * Anchors (mbpriv.h): sid = tid<<1|rev, len = seed length, qpos/tpos = query
 * and target coordinates of the LAST base (query coordinates of a reverse
 * anchor are on the reverse-complemented read). The pair score in the
 * minimap2 paper's form: min(len_i, min(dq, dr)) - (int)(g*dd + 0.5*log2(dd+1))
 * with dd = |dr - dq|; minibwa passes g = 0.8 * 0.01 * min_len = 0.152. */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "mbpriv.h"
#include "minibwa_test.h"

#define BC_PEN 0.152f

static uint64_t bc_state;

static uint64_t bc_rand(void)
{
	uint64_t z = (bc_state += 0x9E3779B97F4A7C15ULL);
	z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
	z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
	return z ^ (z >> 31);
}

static int bc_below(int n) { return (int)(bc_rand() % (uint64_t)n); }

static l2b_ctg_t bc_ctg[2] = { { .name = "c0", .len = 100000, .off = 0 }, { .name = "c1", .len = 50000, .off = 100000 } };
static l2b_t bc_l2b; // only ctg[].off and ctg[].len are read by the code under test

static void bc_init_l2b(void)
{
	memset(&bc_l2b, 0, sizeof(bc_l2b));
	bc_l2b.n_ctg = 2;
	bc_l2b.ctg = bc_ctg;
	bc_l2b.tot_len = 150000;
}

/* an anchor from its start coordinates */
static mb_anchor_t bc_anc(int sid, int qs, int64_t ts, int len)
{
	mb_anchor_t a;
	memset(&a, 0, sizeof(a));
	a.sid = sid, a.len = len, a.qpos = qs + len - 1, a.tpos = ts + len - 1;
	return a;
}

typedef struct {
	int n_u;
	uint64_t *u;
	mb_anchor_t *a; // anchors of chain 0, then chain 1, ...
} bc_chains_t;

static void bc_chain(int n, const mb_anchor_t *in, int max_dist, int bw, int min_sc, bc_chains_t *c)
{
	mb_anchor_t *a = (mb_anchor_t*)malloc((n > 0? n : 1) * sizeof(mb_anchor_t)); // mb_lchain_dp frees its input
	memcpy(a, in, n * sizeof(mb_anchor_t));
	c->u = 0;
	c->a = mb_lchain_dp(0, &bc_l2b, max_dist, max_dist, bw, 25, 5000, min_sc, BC_PEN, n, a, &c->n_u, &c->u);
}

static void bc_free(bc_chains_t *c) { free(c->u); free(c->a); }

/* check chain k: its score and that its anchors are exactly the listed input anchors in order */
static void bc_expect(const bc_chains_t *c, int k, int score, int cnt, const mb_anchor_t *in, const int *idx, const char *what)
{
	int i, off = 0, ok = 1;
	if (!MBT_CHECKF(k < c->n_u, "%s: chain %d missing (%d chains)", what, k, c->n_u)) return;
	for (i = 0; i < k; ++i) off += (int32_t)c->u[i];
	MBT_CHECKF((int)(c->u[k]>>32) == score, "%s: chain %d score %d, expected %d", what, k, (int)(c->u[k]>>32), score);
	if (!MBT_CHECKF((int32_t)c->u[k] == cnt, "%s: chain %d has %d anchors, expected %d", what, k, (int32_t)c->u[k], cnt)) return;
	for (i = 0; i < cnt; ++i) {
		const mb_anchor_t *x = &c->a[off + i], *y = &in[idx[i]];
		if (x->sid != y->sid || x->qpos != y->qpos || x->tpos != y->tpos || x->len != y->len) ok = 0;
	}
	MBT_CHECKF(ok, "%s: chain %d holds the wrong anchors", what, k);
}

/* Hand-worked chains. Each block states the anchors and the arithmetic. */
static int case_lchain_hand(int argc, char *argv[])
{
	bc_chains_t c;
	mb_anchor_t a[8];
	bc_init_l2b();
	{ // three collinear 20-mers, no gaps: 20+20+20
		static const int i0[] = {0, 1, 2};
		a[0] = bc_anc(0, 0, 1000, 20), a[1] = bc_anc(0, 20, 1020, 20), a[2] = bc_anc(0, 40, 1040, 20);
		bc_chain(3, a, 100, 100, 25, &c);
		MBT_CHECK_INT(c.n_u, ==, 1);
		bc_expect(&c, 0, 60, 3, a, i0, "collinear");
		bc_free(&c);
	}
	{ // dq=40, dr=50: dd=10, 20 - (int)(1.52 + 0.5*log2(11)=1.73) = 17 -> 37
		static const int i0[] = {0, 1};
		a[0] = bc_anc(0, 0, 1000, 20), a[1] = bc_anc(0, 40, 1050, 20);
		bc_chain(2, a, 100, 100, 25, &c);
		MBT_CHECK_INT(c.n_u, ==, 1);
		bc_expect(&c, 0, 37, 2, a, i0, "gap 10");
		bc_free(&c);
	}
	{ // a 30 bp deletion between two 75-mers: 75 + 75 - (int)(4.56 + 0.5*log2(31)=2.48) = 143
		static const int i0[] = {0, 1};
		a[0] = bc_anc(0, 0, 5000, 75), a[1] = bc_anc(0, 75, 5105, 75);
		bc_chain(2, a, 150, 100, 25, &c);
		MBT_CHECK_INT(c.n_u, ==, 1);
		bc_expect(&c, 0, 143, 2, a, i0, "deletion 30");
		bc_free(&c);
	}
	{ // overlapping 30-mers 10 bp apart: 30 + min(30, 10) = 40
		static const int i0[] = {0, 1};
		a[0] = bc_anc(0, 0, 1000, 30), a[1] = bc_anc(0, 10, 1010, 30);
		bc_chain(2, a, 100, 100, 25, &c);
		MBT_CHECK_INT(c.n_u, ==, 1);
		bc_expect(&c, 0, 40, 2, a, i0, "overlap");
		bc_free(&c);
	}
	{ // dd=15: joined with bw 20 (20 + 20 - (int)(2.28 + 2) = 36), split with bw 10
		static const int i0[] = {0, 1}, j0[] = {0}, j1[] = {1};
		a[0] = bc_anc(0, 0, 1000, 20), a[1] = bc_anc(0, 40, 1055, 20);
		bc_chain(2, a, 100, 20, 15, &c);
		MBT_CHECK_INT(c.n_u, ==, 1);
		bc_expect(&c, 0, 36, 2, a, i0, "bw 20");
		bc_free(&c);
		bc_chain(2, a, 100, 10, 15, &c);
		MBT_CHECK_INT(c.n_u, ==, 2);
		bc_expect(&c, 0, 20, 1, a, j0, "bw 10 first");
		bc_expect(&c, 1, 20, 1, a, j1, "bw 10 second");
		bc_free(&c);
		bc_chain(2, a, 100, 15, 15, &c); // dd equal to the band is inside it
		MBT_CHECK_INT(c.n_u, ==, 1);
		bc_expect(&c, 0, 36, 2, a, i0, "bw 15");
		bc_free(&c);
	}
	{ // a step that costs more than it gains: A (60) -> B (20, dd 150, dg 40) scores 20 - (int)(22.8 + 3.6) = -6,
	  // then B -> C (60, collinear) +60. The best chain runs through B: 60 - 6 + 60 = 114 (A -> C directly:
	  // 60 + 60 - 26 = 94) and is reported whole, not cut at the dip
		static const int i0[] = {0, 1, 2};
		a[0] = bc_anc(0, 0, 1000, 60), a[1] = bc_anc(0, 80, 1230, 20), a[2] = bc_anc(0, 100, 1250, 60);
		bc_chain(3, a, 300, 200, 25, &c);
		MBT_CHECK_INT(c.n_u, ==, 1);
		bc_expect(&c, 0, 114, 3, a, i0, "negative step");
		bc_free(&c);
	}
	{ // 150 bp apart on both sequences: joined when max_dist (and so the window) is 200, not at 100 (bw 100)
		static const int i0[] = {0, 1}, j0[] = {0}, j1[] = {1};
		a[0] = bc_anc(0, 0, 1000, 20), a[1] = bc_anc(0, 150, 1150, 20);
		bc_chain(2, a, 200, 100, 15, &c);
		MBT_CHECK_INT(c.n_u, ==, 1);
		bc_expect(&c, 0, 40, 2, a, i0, "max_dist 200");
		bc_free(&c);
		bc_chain(2, a, 100, 100, 15, &c);
		MBT_CHECK_INT(c.n_u, ==, 2);
		bc_expect(&c, 0, 20, 1, a, j0, "max_dist 100 first");
		bc_expect(&c, 1, 20, 1, a, j1, "max_dist 100 second");
		bc_free(&c);
	}
	{ // min_sc: a 40-point chain is reported at 40 and dropped at 41
		static const int i0[] = {0, 1};
		a[0] = bc_anc(0, 0, 1000, 20), a[1] = bc_anc(0, 20, 1020, 20);
		bc_chain(2, a, 100, 100, 40, &c);
		bc_expect(&c, 0, 40, 2, a, i0, "min_sc 40");
		bc_free(&c);
		bc_chain(2, a, 100, 100, 41, &c);
		MBT_CHECK(c.n_u == 0 && c.a == 0);
		bc_free(&c);
	}
	{ // query order against target order: anchors that go backwards on the target do not chain
		a[0] = bc_anc(0, 40, 990, 20), a[1] = bc_anc(0, 0, 1000, 20); // sorted by target
		bc_chain(2, a, 100, 100, 15, &c);
		MBT_CHECK_INT(c.n_u, ==, 2);
		bc_free(&c);
	}
	{ // strands and contigs never mix; chains come back ordered by contig, strand, position
		static const int f0[] = {0, 1, 2}, r0[] = {3, 4}, c1[] = {5, 6};
		a[0] = bc_anc(0, 0, 5000, 20), a[1] = bc_anc(0, 20, 5020, 20), a[2] = bc_anc(0, 40, 5040, 20);
		a[3] = bc_anc(1, 10, 3000, 20), a[4] = bc_anc(1, 30, 3020, 20);
		a[5] = bc_anc(2, 60, 100, 25), a[6] = bc_anc(2, 85, 125, 25);
		bc_chain(7, a, 100, 100, 25, &c);
		MBT_CHECK_INT(c.n_u, ==, 3);
		bc_expect(&c, 0, 60, 3, a, f0, "contig 0 forward");
		bc_expect(&c, 1, 40, 2, a, r0, "contig 0 reverse");
		bc_expect(&c, 2, 50, 2, a, c1, "contig 1 forward");
		bc_free(&c);
	}
	{ // no anchors
		mb_anchor_t *p = 0;
		int n_u = -1;
		uint64_t *u = (uint64_t*)1;
		p = mb_lchain_dp(0, &bc_l2b, 100, 100, 100, 25, 5000, 25, BC_PEN, 0, 0, &n_u, &u);
		MBT_CHECK(p == 0 && n_u == 0 && u == 0);
	}
	return 0;
}

/* the documented pair score, written independently of lchain.c */
static int bc_pair_sc(const mb_anchor_t *aj, const mb_anchor_t *ai, int max_dist, int bw)
{
	int64_t dq = ai->qpos - aj->qpos, dr = ai->tpos - aj->tpos, dd, dg, sc;
	if (ai->sid != aj->sid || dq <= 0 || dr <= 0) return INT32_MIN;
	if (dq > max_dist + ai->len || dr >= max_dist + ai->len) return INT32_MIN;
	dd = dr > dq? dr - dq : dq - dr;
	if (dd > bw) return INT32_MIN;
	dg = dr < dq? dr : dq;
	sc = ai->len < dg? ai->len : dg;
	if (dd) sc -= (int)(BC_PEN * (float)dd + .5f * mb_log2(dd + 1));
	return (int)sc;
}

static int bc_cmp_anchor(const void *x, const void *y)
{
	const mb_anchor_t *a = (const mb_anchor_t*)x, *b = (const mb_anchor_t*)y;
	if (a->sid != b->sid) return a->sid < b->sid? -1 : 1;
	if (a->tpos != b->tpos) return a->tpos < b->tpos? -1 : 1;
	return a->qpos < b->qpos? -1 : a->qpos > b->qpos? 1 : 0;
}

/* Random anchors: one true chain of 10-60 seeds along a slowly drifting diagonal plus
 * noise seeds. The best chain mb_lchain_dp reports scores the brute-force optimum,
 * its anchors form a valid chain whose pairwise scores add up to that, and no
 * anchor is used twice. args: seed n_iter */
static int case_lchain_random(int argc, char *argv[])
{
	int it, n_iter;
	bc_init_l2b();
	if (argc < 2) return MBT_FAIL("usage: lchain_random <seed> <n_iter>"), 1;
	bc_state = strtoull(argv[0], 0, 10);
	n_iter = atoi(argv[1]);
	for (it = 0; it < n_iter; ++it) {
		int n_true = 10 + bc_below(50), n_noise = bc_below(80), n = 0, i, j, k, best = 0, top = -1, off, top_off = 0;
		int qs = bc_below(50), sid = bc_below(4), max_dist = 100 + bc_below(200), bw = 50 + bc_below(100);
		int64_t ts = 1000 + bc_below(20000);
		mb_anchor_t *a = (mb_anchor_t*)malloc((n_true + n_noise) * sizeof(mb_anchor_t));
		int *f = (int*)malloc((n_true + n_noise) * sizeof(int));
		bc_chains_t c;
		for (i = 0; i < n_true; ++i) {
			int len = 19 + bc_below(30), gap = 1 + bc_below(40), d = bc_below(10) == 0? bc_below(7) - 3 : 0;
			a[n++] = bc_anc(sid, qs, ts, len);
			qs += len + gap, ts += len + gap + d;
			if (ts < 0) ts = 0;
		}
		for (i = 0; i < n_noise; ++i) {
			int s = bc_below(4);
			a[n++] = bc_anc(s, bc_below(qs + 1), bc_below(s < 2? 90000 : 40000), 19 + bc_below(20));
		}
		qsort(a, n, sizeof(mb_anchor_t), bc_cmp_anchor);
		if (max_dist < bw) max_dist = bw; // mb_lchain_dp never looks less far than the band
		for (i = 0; i < n; ++i) { // brute force: f[i] = max(len_i, max_j f[j] + sc(j, i))
			f[i] = a[i].len;
			for (j = 0; j < i; ++j) {
				int sc = bc_pair_sc(&a[j], &a[i], max_dist, bw);
				if (sc != INT32_MIN && f[j] + sc > f[i]) f[i] = f[j] + sc;
			}
			if (f[i] > best) best = f[i];
		}
		bc_chain(n, a, max_dist, bw, 25, &c);
		for (k = 0, off = 0; k < c.n_u; ++k) {
			if ((int)(c.u[k]>>32) > (top < 0? -1 : (int)(c.u[top]>>32))) top = k, top_off = off;
			for (i = 1; i < (int32_t)c.u[k]; ++i) {
				const mb_anchor_t *p = &c.a[off + i - 1], *q = &c.a[off + i];
				if (p->sid != q->sid || p->qpos >= q->qpos || p->tpos >= q->tpos) {
					MBT_FAIL("iter %d: chain %d is not collinear at anchor %d", it, k, i);
					break;
				}
			}
			off += (int32_t)c.u[k];
		}
		for (i = 0; i < off; ++i) // each input anchor in at most one chain
			for (j = i + 1; j < off; ++j)
				if (c.a[i].sid == c.a[j].sid && c.a[i].qpos == c.a[j].qpos && c.a[i].tpos == c.a[j].tpos && c.a[i].len == c.a[j].len)
					MBT_FAIL("iter %d: an anchor appears twice", it), i = j = off;
		if (MBT_CHECKF(top >= 0, "iter %d: no chain", it)) {
			int s = c.a[top_off].len;
			MBT_CHECKF((int)(c.u[top]>>32) == best, "iter %d: best chain %d, brute force %d", it, (int)(c.u[top]>>32), best);
			for (i = 1; i < (int32_t)c.u[top]; ++i)
				s += bc_pair_sc(&c.a[top_off + i - 1], &c.a[top_off + i], max_dist, bw);
			MBT_CHECKF(s == best, "iter %d: the best chain's anchors score %d, reported %d", it, s, best);
		}
		bc_free(&c);
		free(f); free(a);
	}
	return 0;
}

/* mb_cal_mblen: mlen counts bases covered by seeds (overlaps once), blen spans the
 * longer of the two gaps. */
static int case_cal_mblen(int argc, char *argv[])
{
	mb_anchor_t a[4];
	int32_t blen, mlen;
	a[0] = bc_anc(0, 0, 1000, 20), a[1] = bc_anc(0, 40, 1050, 20); // gap: tl 50, ql 40
	mlen = mb_cal_mblen(2, a, &blen);
	MBT_CHECK(mlen == 40 && blen == 70);
	a[0] = bc_anc(0, 0, 1000, 30), a[1] = bc_anc(0, 10, 1010, 30); // overlap by 20
	mlen = mb_cal_mblen(2, a, &blen);
	MBT_CHECK(mlen == 40 && blen == 40);
	a[0] = bc_anc(0, 0, 1000, 75), a[1] = bc_anc(0, 76, 1076, 74); // a mismatch at 75
	mlen = mb_cal_mblen(2, a, &blen);
	MBT_CHECK(mlen == 149 && blen == 150);
	a[0] = bc_anc(0, 0, 5000, 75), a[1] = bc_anc(0, 75, 5105, 75); // a 30 bp deletion
	mlen = mb_cal_mblen(2, a, &blen);
	MBT_CHECK(mlen == 150 && blen == 180);
	mlen = mb_cal_mblen(1, a, &blen);
	MBT_CHECK(mlen == 75 && blen == 75);
	mlen = mb_cal_mblen(0, a, &blen);
	MBT_CHECK(mlen == 0 && blen == 0);
	return 0;
}

/* mb_gen_hit: chains become hits sorted by score, without a base alignment; coordinates
 * come from the first and last anchors, query coordinates of reverse hits are on the read's
 * forward strand. (The internal id, parent and anchor-offset fields are not checked.) */
static int case_gen_hit(int argc, char *argv[])
{
	bc_chains_t c;
	mb_anchor_t a[5];
	mb_hit_t *h;
	int i, qlen = 200;
	bc_init_l2b();
	a[0] = bc_anc(1, 10, 3000, 20), a[1] = bc_anc(1, 30, 3020, 20); // reverse strand, score 40
	a[2] = bc_anc(0, 0, 5000, 20), a[3] = bc_anc(0, 20, 5020, 20), a[4] = bc_anc(0, 40, 5040, 20); // forward, 60
	qsort(a, 5, sizeof(mb_anchor_t), bc_cmp_anchor);
	bc_chain(5, a, 100, 100, 25, &c);
	MBT_REQUIRE(c.n_u == 2);
	h = mb_gen_hit(0, 12345, qlen, &bc_l2b, c.n_u, c.u, c.a);
	MBT_REQUIRE(h != 0);
	MBT_CHECK(h[0].score == 60 && h[0].score0 == 60 && h[0].cnt == 3 && h[1].score == 40 && h[1].cnt == 2);
	MBT_CHECK(h[0].tid == 0 && h[0].rev == 0 && h[0].ts == 5000 && h[0].te == 5060 && h[0].qs == 0 && h[0].qe == 60);
	MBT_CHECK(h[1].tid == 0 && h[1].rev == 1 && h[1].ts == 3000 && h[1].te == 3040 && h[1].qs == qlen - 50 && h[1].qe == qlen - 10);
	MBT_CHECK(h[0].mlen == 60 && h[0].blen == 60 && h[1].mlen == 40 && h[1].blen == 40);
	for (i = 0; i < 2; ++i)
		MBT_CHECK(h[i].p == 0);
	free(h);
	MBT_CHECK(mb_gen_hit(0, 1, qlen, &bc_l2b, 0, c.u, c.a) == 0);
	bc_free(&c);
	return 0;
}

static void bc_hit_from(mb_hit_t *h, int as, int cnt, int score, int qlen, const mb_anchor_t *a)
{
	memset(h, 0, sizeof(*h));
	h->as = as, h->cnt = cnt, h->score = h->score0 = score;
	h->tid = a[as].sid >> 1, h->rev = a[as].sid & 1;
	h->ts = a[as].tpos + 1 - a[as].len, h->te = a[as + cnt - 1].tpos + 1;
	h->qs = h->rev? qlen - (a[as + cnt - 1].qpos + 1) : a[as].qpos + 1 - a[as].len;
	h->qe = h->rev? qlen - (a[as].qpos + 1 - a[as].len) : a[as + cnt - 1].qpos + 1;
}

/* mb_split_hit: the first n anchors stay, the rest form a new hit with a proportional
 * share of the score; the piece split off a primary becomes a primary of its own once the
 * hits are synced (mb_sync_hits), and is not the SAM primary. */
static int case_split_hit(int argc, char *argv[])
{
	mb_anchor_t a[4];
	mb_hit_t r, r2, hs[2];
	bc_init_l2b();
	a[0] = bc_anc(0, 0, 1000, 20), a[1] = bc_anc(0, 20, 1020, 20), a[2] = bc_anc(0, 50, 1060, 20), a[3] = bc_anc(0, 70, 1080, 20);
	bc_hit_from(&r, 0, 4, 77, 100, a);
	r.id = r.parent = 0, r.sam_pri = 1;
	mb_split_hit(&r, &r2, 1, 100, a, &bc_l2b);
	MBT_CHECK(r.cnt == 1 && r2.cnt == 3);
	MBT_CHECK_INT(r2.score, ==, 58); // 77 * 3/4 = 57.75, rounded
	MBT_CHECK_INT(r.score, ==, 19);
	MBT_CHECK(r.ts == 1000 && r.te == 1020 && r.qs == 0 && r.qe == 20);
	MBT_CHECK(r2.ts == 1020 && r2.te == 1100 && r2.qs == 20 && r2.qe == 90);
	MBT_CHECK(r2.sam_pri == 0 && r2.p == 0);
	MBT_CHECK(r.split == 1 && r2.split == 2);
	MBT_CHECK(r.mlen == 20 && r2.mlen == 60 && r2.blen == 80);
	hs[0] = r, hs[1] = r2;
	mb_sync_hits(0, 2, hs);
	MBT_CHECK(hs[0].parent == hs[0].id && hs[1].parent == hs[1].id && hs[0].id != hs[1].id);
	MBT_CHECK(hs[0].sam_pri == 1 && hs[1].sam_pri == 0);
	bc_hit_from(&r, 0, 4, 80, 100, a); // a secondary keeps its parent; out-of-range n does nothing
	r.id = 1, r.parent = 0;
	mb_split_hit(&r, &r2, 2, 100, a, &bc_l2b);
	MBT_CHECK(r2.parent == 0 && r2.score == 40 && r.score == 40);
	bc_hit_from(&r, 0, 4, 80, 100, a);
	r2.cnt = -5;
	mb_split_hit(&r, &r2, 4, 100, a, &bc_l2b);
	mb_split_hit(&r, &r2, 0, 100, a, &bc_l2b);
	MBT_CHECK(r.cnt == 4 && r.score == 80 && r2.cnt == -5);
	{ // reverse strand: query coordinates on the read's forward strand
		mb_anchor_t b[2];
		b[0] = bc_anc(1, 10, 3000, 20), b[1] = bc_anc(1, 30, 3020, 20);
		bc_hit_from(&r, 0, 2, 40, 200, b);
		mb_split_hit(&r, &r2, 1, 200, b, &bc_l2b);
		MBT_CHECK(r.qs == 170 && r.qe == 190 && r2.qs == 150 && r2.qe == 170 && r2.rev == 1);
	}
	return 0;
}

/* mb_squeeze_a packs the anchors of the remaining hits into the first n slots; each hit
 * still points at its own anchors, in order (where in [0, n) they go is not checked) */
static int case_squeeze_a(int argc, char *argv[])
{
	mb_anchor_t a[8];
	mb_hit_t h[2];
	int i, n;
	for (i = 0; i < 8; ++i) a[i] = bc_anc(0, i * 10, 1000 + i * 10, 10);
	memset(h, 0, sizeof(h));
	h[0].as = 5, h[0].cnt = 2; // anchors 5, 6 (qpos 59, 69)
	h[1].as = 1, h[1].cnt = 2; // anchors 1, 2 (qpos 19, 29)
	n = mb_squeeze_a(0, 2, h, a);
	MBT_CHECK_INT(n, ==, 4);
	MBT_REQUIRE(h[0].as >= 0 && h[0].as + 2 <= n && h[1].as >= 0 && h[1].as + 2 <= n);
	MBT_CHECK(h[0].as + 2 <= h[1].as || h[1].as + 2 <= h[0].as);
	MBT_CHECK(a[h[0].as].qpos == 59 && a[h[0].as + 1].qpos == 69);
	MBT_CHECK(a[h[1].as].qpos == 19 && a[h[1].as + 1].qpos == 29);
	return 0;
}

/* mb_sync_hits: ids become array indices, parents follow their hit to its new index,
 * a provisional primary becomes its own parent, a parent that is gone becomes unset */
static int case_sync_hits(int argc, char *argv[])
{
	mb_hit_t h[4];
	memset(h, 0, sizeof(h));
	h[0].id = 3, h[0].parent = 3, h[0].qs = 50;  // a primary that was hit 3
	h[1].id = 0, h[1].parent = 3;                // its secondary
	h[2].id = -1, h[2].parent = MB_PARENT_TMP_PRI, h[2].qs = 10; // split off a primary
	h[3].id = 5, h[3].parent = 4;                // hit 4 no longer exists
	mb_sync_hits(0, 4, h);
	MBT_CHECK(h[0].id == 0 && h[1].id == 1 && h[2].id == 2 && h[3].id == 3);
	MBT_CHECK(h[0].parent == 0 && h[1].parent == 0 && h[2].parent == 2 && h[3].parent == MB_PARENT_UNSET);
	MBT_CHECK(h[0].sam_pri == 1 && h[1].sam_pri == 0 && h[2].sam_pri == 0);
	return 0;
}

static const mbt_case_t cases[] = {
	{ "lchain_hand", case_lchain_hand },
	{ "lchain_random", case_lchain_random },
	{ "cal_mblen", case_cal_mblen },
	{ "gen_hit", case_gen_hit },
	{ "split_hit", case_split_hit },
	{ "squeeze_a", case_squeeze_a },
	{ "sync_hits", case_sync_hits },
	{ 0, 0 }
};

MBT_MAIN(cases)
