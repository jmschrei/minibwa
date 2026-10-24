/* Area B driver: hit processing and mapping quality in map-algo.c
 * (mb_set_parent, mb_set_sam_pri, mb_select_sub, mb_hit_sort,
 * mb_filter_hits, mb_set_mapq, mb_cal_high_cov), the presets and the
 * adaptive parameters in options.c, the inline helpers in mbpriv.h and the
 * radix sorts. Expected values are worked by hand from the formulas named in
 * each case; the MAPQ arithmetic is spelled out next to each expectation. */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include "mbpriv.h"
#include "minibwa_test.h"

// defined in map-algo.c without a prototype in a header
int32_t mb_cal_high_cov(void *km, int32_t n, const mb_sai_t *sai, int32_t max_occ);

static mb_extra_t *bh_extra(int dp_max, int dp_max2)
{
	mb_extra_t *p = (mb_extra_t*)calloc(1, sizeof(mb_extra_t) + 16);
	p->dp_max = p->dp_max0 = dp_max, p->dp_max2 = dp_max2;
	return p;
}

static void bh_hit(mb_hit_t *h, int id, int qs, int qe, int score)
{
	memset(h, 0, sizeof(*h));
	h->id = id, h->parent = MB_PARENT_UNSET;
	h->qs = qs, h->qe = qe, h->score = h->score0 = score;
	h->ts = 1000 + id * 1000, h->te = h->ts + (qe - qs);
	h->cnt = 1;
}

/* mb_set_parent: a hit becomes a secondary of an earlier primary when the query overlap
 * over the shorter hit, minus its uncovered fraction, exceeds mask_level (strictly). */
static int case_set_parent(int argc, char *argv[])
{
	mb_hit_t h[6];
	int i;
	bh_hit(&h[0], 0, 0, 100, 100);
	bh_hit(&h[1], 1, 10, 100, 80);  // inside h0: secondary of 0
	bh_hit(&h[2], 2, 100, 150, 50); // no overlap: primary
	bh_hit(&h[3], 3, 60, 140, 40);  // 40/80 = 0.5 with h0 (not > 0.5); 40/50 with h2: secondary of 2
	bh_hit(&h[4], 4, 0, 200, 30);   // 100/100 - 50 uncovered/200 = 0.75 with h0: secondary of 0
	bh_hit(&h[5], 5, 120, 300, 20); // 30/50 - 150/180 < 0.5: primary
	mb_set_parent(0, 0.5f, INT_MAX, 6, h, 14, 0);
	MBT_CHECK(h[0].parent == 0 && h[1].parent == 0 && h[2].parent == 2 && h[3].parent == 2 && h[4].parent == 0 && h[5].parent == 5);
	MBT_CHECK(h[0].subsc == 80 && h[2].subsc == 40 && h[5].subsc == 0);
	for (i = 0; i < 6; ++i) MBT_CHECK_INT(h[i].id, ==, i);
	// mask_len caps the uncovered length a secondary may have: h4 leaves 50 bp uncovered
	for (i = 0; i < 6; ++i) h[i].parent = MB_PARENT_UNSET, h[i].subsc = 0;
	mb_set_parent(0, 0.5f, 40, 6, h, 14, 0);
	MBT_CHECK(h[4].parent == 4 && h[1].parent == 0);
	// a higher mask level keeps h1 (overlap fraction 1.0) but frees h3 (0.8 at best)
	for (i = 0; i < 6; ++i) h[i].parent = MB_PARENT_UNSET, h[i].subsc = 0;
	mb_set_parent(0, 0.85f, INT_MAX, 6, h, 14, 0);
	MBT_CHECK(h[1].parent == 0 && h[3].parent == 3);
	// dp_max2 and n_sub from base-aligned secondaries; a hit identical to the parent does not count
	{
		mb_hit_t g[4];
		bh_hit(&g[0], 0, 0, 150, 150), g[0].p = bh_extra(300, 0);
		bh_hit(&g[1], 1, 0, 150, 140), g[1].p = bh_extra(290, 0); // within sub_diff 14
		bh_hit(&g[2], 2, 0, 150, 120), g[2].p = bh_extra(250, 0); // 50 below
		g[3] = g[0], g[3].id = 3, g[3].p = bh_extra(299, 0); // same coordinates as g0
		mb_set_parent(0, 0.5f, INT_MAX, 4, g, 14, 0);
		MBT_CHECK(g[1].parent == 0 && g[2].parent == 0 && g[3].parent == 0);
		MBT_CHECK_INT(g[0].p->dp_max2, ==, 290);
		MBT_CHECK_INT(g[0].n_sub, ==, 1);
		MBT_CHECK_INT(g[0].subsc, ==, 150);
		for (i = 0; i < 4; ++i) free(g[i].p);
	}
	return 0;
}

/* mb_set_sam_pri: the first primary, or the primary with the smallest query start */
static int case_set_sam_pri(int argc, char *argv[])
{
	mb_hit_t h[4];
	int k;
	bh_hit(&h[0], 0, 50, 150, 100), h[0].parent = 0;
	bh_hit(&h[1], 1, 50, 150, 90), h[1].parent = 0;
	bh_hit(&h[2], 2, 10, 50, 40), h[2].parent = 2;
	bh_hit(&h[3], 3, 0, 40, 30), h[3].parent = 2; // smallest qs, but a secondary
	k = mb_set_sam_pri(4, h, 0);
	MBT_CHECK(k == 0 && h[0].sam_pri && !h[1].sam_pri && !h[2].sam_pri && !h[3].sam_pri);
	k = mb_set_sam_pri(4, h, 1);
	MBT_CHECK(k == 2 && !h[0].sam_pri && h[2].sam_pri && !h[3].sam_pri);
	MBT_CHECK_INT(mb_set_sam_pri(0, h, 0), ==, -1);
	return 0;
}

/* mb_select_sub: a secondary survives if its chain score is at least pri_ratio times its
 * parent's or within min_diff of it, it is not a copy of its parent, and fewer than
 * best_n secondaries have survived; survivors are renumbered */
static int case_select_sub(int argc, char *argv[])
{
	mb_hit_t h[7];
	int n;
	bh_hit(&h[0], 0, 0, 100, 100), h[0].parent = 0;
	bh_hit(&h[1], 1, 0, 100, 60), h[1].parent = 0;  // 60 >= 50: kept
	bh_hit(&h[2], 2, 0, 100, 40), h[2].parent = 0;  // 40 < 50 and 40 + 38 < 100: dropped
	h[3] = h[0], h[3].id = 3, h[3].score = 70;      // a copy of h0: dropped
	bh_hit(&h[4], 4, 100, 150, 50), h[4].parent = 4;
	bh_hit(&h[5], 5, 100, 150, 30), h[5].parent = 4; // 30 >= 25: kept
	bh_hit(&h[6], 6, 0, 100, 99), h[6].parent = 0;   // over the best_n cap of 2
	n = 7;
	mb_select_sub(0, 0.5f, 38, 2, &n, h);
	MBT_REQUIRE(n == 4);
	MBT_CHECK(h[0].score == 100 && h[1].score == 60 && h[2].score == 50 && h[3].score == 30);
	MBT_CHECK(h[0].id == 0 && h[1].id == 1 && h[2].id == 2 && h[3].id == 3);
	MBT_CHECK(h[0].parent == 0 && h[1].parent == 0 && h[2].parent == 2 && h[3].parent == 2);
	// min_diff rescues a secondary below the ratio: 60 + 38 < 100 is dropped at ratio 0.9, 99 kept
	bh_hit(&h[0], 0, 0, 100, 100), h[0].parent = 0;
	bh_hit(&h[1], 1, 0, 100, 60), h[1].parent = 0;
	bh_hit(&h[2], 2, 0, 100, 63), h[2].parent = 0; // 63 + 38 >= 100: kept
	n = 3;
	mb_select_sub(0, 0.9f, 38, 50, &n, h);
	MBT_CHECK(n == 2 && h[1].score == 63 && h[1].parent == 0);
	// an inversion hit is always kept; ratio 0 keeps everything
	bh_hit(&h[0], 0, 0, 100, 100), h[0].parent = 0;
	bh_hit(&h[1], 1, 0, 100, 1), h[1].parent = 0, h[1].inv = 1;
	bh_hit(&h[2], 2, 0, 100, 2), h[2].parent = 0;
	n = 3;
	mb_select_sub(0, 0.5f, 0, 50, &n, h);
	MBT_CHECK(n == 2 && h[1].inv);
	bh_hit(&h[2], 2, 0, 100, 2), h[2].parent = 0;
	n = 3;
	mb_select_sub(0, 0.0f, 0, 0, &n, h);
	MBT_CHECK(n == 3);
	return 0;
}

/* mb_hit_sort: by DP score (chain score without a base alignment), highest first; hits
 * with a negative anchor count are dropped unless they are inversions */
static int case_hit_sort(int argc, char *argv[])
{
	mb_hit_t h[5];
	int n = 5;
	bh_hit(&h[0], 0, 0, 10, 1), h[0].p = bh_extra(100, 0), h[0].hash = 7;
	bh_hit(&h[1], 1, 0, 10, 1), h[1].p = bh_extra(300, 0);
	bh_hit(&h[2], 2, 0, 10, 200);
	bh_hit(&h[3], 3, 0, 10, 1), h[3].p = bh_extra(500, 0), h[3].cnt = -1;
	bh_hit(&h[4], 4, 0, 10, 1), h[4].p = bh_extra(50, 0), h[4].cnt = -1, h[4].inv = 1;
	mb_hit_sort(0, &n, h);
	MBT_REQUIRE(n == 4);
	MBT_CHECK(h[0].p && h[0].p->dp_max == 300);
	MBT_CHECK(h[1].p == 0 && h[1].score == 200);
	MBT_CHECK(h[2].p && h[2].p->dp_max == 100 && h[2].hash == 7);
	MBT_CHECK(h[3].inv && h[3].p->dp_max == 50);
	free(h[0].p); free(h[2].p); free(h[3].p);
	return 0;
}

/* mb_filter_hits: an aligned hit needs mlen >= -m and DP max >= -s times -A */
static int case_filter_hits(int argc, char *argv[])
{
	mb_opt_t opt;
	mb_hit_t h[5];
	int n = 5;
	mb_opt_init(&opt);
	MBT_REQUIRE(opt.min_chain_score == 25 && opt.min_dp_max == 30 && opt.a == 2);
	bh_hit(&h[0], 0, 0, 10, 1), h[0].mlen = 24, h[0].p = bh_extra(100, 0);
	bh_hit(&h[1], 1, 0, 10, 1), h[1].mlen = 25, h[1].p = bh_extra(59, 0);
	bh_hit(&h[2], 2, 0, 10, 1), h[2].mlen = 25, h[2].p = bh_extra(60, 0);
	bh_hit(&h[3], 3, 0, 10, 1), h[3].mlen = 0; // no base alignment: kept
	bh_hit(&h[4], 4, 0, 10, 1), h[4].flt = 1;  // flagged: dropped
	mb_filter_hits(&opt, 150, &n, h);
	MBT_REQUIRE(n == 2);
	MBT_CHECK(h[0].id == 2 && h[1].id == 3);
	free(h[0].p);
	return 0;
}

static void bh_mq_hit(mb_hit_t *h, int score, int score0, int subsc, int mlen, int blen, int dp_max, int dp_max2, int n_sub, int seed_ratio)
{
	bh_hit(h, 0, 0, 100, score);
	h->parent = 0, h->score0 = score0, h->subsc = subsc, h->mlen = mlen, h->blen = blen;
	h->n_sub = n_sub, h->seed_ratio = seed_ratio;
	if (dp_max >= 0) h->p = bh_extra(dp_max, dp_max2);
}

static int bh_mapq(int qlen, int is_sr, int max_sr_len, int min_chain_sc, mb_hit_t *h)
{
	int q;
	mb_set_mapq(0, qlen, 1, h, min_chain_sc, 2, is_sr, max_sr_len);
	q = h->mapq;
	free(h->p);
	h->p = 0;
	return q;
}

/* mb_set_mapq with hand-computed values (a = 2):
 * short-read (bwa-mem): x = (blen < 50? 1 : 3/ln(blen)) * identity^2,
 *   mapq_sr = (int)(6.02 * x^2 * (dp_max - dp_max2) / a + .499)
 * long-read (minimap2): x = dp_max2/dp_max (times subsc/score0 if subsc > score0),
 *   mapq_lr = (int)(pen * identity * 40 * (1 - x^2) * ln(dp_max/a))
 * blend for qlen >= -l: mapq_lr - (mapq_lr - mapq_sr) * 2^(1 - qlen/l), rounded
 * without a second hit: (int)(pen * identity * 40 * (1 - max(subsc,-m)/score0) * ln(dp_max/a))
 * then minus (int)(4.343 ln(n_sub+1) + .499), floor 0, times (seed_ratio/50)^2 below 50,
 * cap 60, and at least 1 if dp_max > dp_max2. */
static int case_set_mapq(int argc, char *argv[])
{
	mb_hit_t h, hs[3];
	// sr: x = 3/ln150 = 0.5987, 6.02 x^2 10/2 = 10.79 -> 11; n_sub 1: -3 -> 8
	bh_mq_hit(&h, 150, 150, 75, 150, 150, 300, 290, 1, 255);
	MBT_CHECK_INT(bh_mapq(150, 1, 325, 25, &h), ==, 8);
	// sr, 100 bp: 3/ln100 = 0.6514 -> 12.77 -> 13, -3 -> 10
	bh_mq_hit(&h, 100, 100, 50, 100, 100, 200, 190, 1, 255);
	MBT_CHECK_INT(bh_mapq(100, 1, 325, 25, &h), ==, 10);
	// sr with mismatches: identity 140/150; x = 0.5987 * 0.8711 = 0.5216; 6.02 x^2 20/2 = 16.4 -> 16
	bh_mq_hit(&h, 150, 150, 75, 140, 150, 270, 250, 0, 255);
	MBT_CHECK_INT(bh_mapq(150, 1, 325, 25, &h), ==, 16);
	// lr (-l < 0): identity 0.95, x = 1700/1800; 0.95*40*(1-0.892)*ln900 = 27.9 -> 27
	bh_mq_hit(&h, 900, 900, 300, 950, 1000, 1800, 1700, 0, 255);
	MBT_CHECK_INT(bh_mapq(1000, 0, -1, 40, &h), ==, 27);
	// lr with subsc > score0: x = 0.944 * 1000/900 > 1, mapq < 0 -> 0 -> raised to 1
	bh_mq_hit(&h, 900, 900, 1000, 950, 1000, 1800, 1700, 0, 255);
	MBT_CHECK_INT(bh_mapq(1000, 0, -1, 40, &h), ==, 1);
	// blend at 400 bp: sr 39 (3/ln400, diff 52), lr 30 (x = 748/800, ln400); 30 + 9 * 2^(1-400/325) = 37.7 -> 38
	bh_mq_hit(&h, 400, 400, 200, 400, 400, 800, 748, 0, 255);
	MBT_CHECK_INT(bh_mapq(400, 0, 325, 25, &h), ==, 38);
	// blend at 500 bp: sr 36, lr 25 -> 25 + 11 * 2^(1-500/325) = 32.6 -> 33
	bh_mq_hit(&h, 500, 500, 250, 500, 500, 1000, 948, 0, 255);
	MBT_CHECK_INT(bh_mapq(500, 0, 325, 25, &h), ==, 33);
	// long-read mode below -l uses the sr formula: 3/ln200 -> 10, n_sub 1: -3 -> 7
	bh_mq_hit(&h, 200, 200, 100, 200, 200, 400, 390, 1, 255);
	MBT_CHECK_INT(bh_mapq(200, 0, 325, 25, &h), ==, 7);
	// no second hit: subsc raised to -m 25; 40 * (1 - 25/31) * ln31 = 26.6 -> 26
	bh_mq_hit(&h, 31, 31, 0, 31, 31, 62, 0, 0, 255);
	MBT_CHECK_INT(bh_mapq(150, 1, 325, 25, &h), ==, 26);
	// few matching bases: pen = 10*50/1000 = 0.5; 0.5 * 40 * 0.5 * ln50 = 39.1 -> 39
	bh_mq_hit(&h, 50, 50, 0, 50, 50, 100, 0, 0, 255);
	MBT_CHECK_INT(bh_mapq(1000, 0, 325, 25, &h), ==, 39);
	// no base alignment: 40 * (1 - 75/150) * ln150 = 100.2 -> cap 60; subsc 140: 13.4 -> 13
	bh_mq_hit(&h, 150, 150, 75, 0, 0, -1, 0, 0, 255);
	MBT_CHECK_INT(bh_mapq(150, 1, 325, 25, &h), ==, 60);
	bh_mq_hit(&h, 150, 150, 140, 0, 0, -1, 0, 0, 255);
	MBT_CHECK_INT(bh_mapq(150, 1, 325, 25, &h), ==, 13);
	// seed ratio 25 scales by (25/50)^2 before the cap: 100 * 0.25 = 25
	bh_mq_hit(&h, 150, 150, 75, 0, 0, -1, 0, 0, 25);
	MBT_CHECK_INT(bh_mapq(150, 1, 325, 25, &h), ==, 25);
	// dp_max one above dp_max2 with many suboptimal hits: 0, raised to 1; equal scores stay 0
	bh_mq_hit(&h, 150, 150, 75, 150, 150, 300, 299, 10, 255);
	MBT_CHECK_INT(bh_mapq(150, 1, 325, 25, &h), ==, 1);
	bh_mq_hit(&h, 150, 150, 75, 150, 150, 300, 300, 0, 255);
	MBT_CHECK_INT(bh_mapq(150, 1, 325, 25, &h), ==, 0);
	// secondaries get 0
	bh_mq_hit(&h, 150, 150, 75, 150, 150, 300, 200, 0, 255);
	h.id = 1;
	MBT_CHECK_INT(bh_mapq(150, 1, 325, 25, &h), ==, 0);
	// an inversion between two primaries takes the smaller of their MAPQs
	bh_mq_hit(&hs[0], 150, 150, 75, 0, 0, -1, 0, 0, 255), hs[0].ts = 1000;
	bh_mq_hit(&hs[1], 150, 150, 140, 0, 0, -1, 0, 0, 255), hs[1].ts = 3000, hs[1].id = hs[1].parent = 1;
	bh_mq_hit(&hs[2], 10, 10, 0, 0, 0, -1, 0, 0, 255), hs[2].ts = 2000, hs[2].id = hs[2].parent = 2, hs[2].inv = 1;
	mb_set_mapq(0, 150, 3, hs, 25, 2, 1, 325);
	MBT_CHECK(hs[0].mapq == 60 && hs[1].mapq == 13 && hs[2].mapq == 13);
	return 0;
}

/* mb_cal_high_cov: query bases covered by seeds occurring more than max_occ times */
static int case_high_cov(int argc, char *argv[])
{
	mb_sai_t s[5];
	memset(s, 0, sizeof(s));
	s[0].info = (uint64_t)0<<32 | 30, s[0].size = 300;
	s[1].info = (uint64_t)20<<32 | 50, s[1].size = 500;
	s[2].info = (uint64_t)70<<32 | 75, s[2].size = 1000;
	s[3].info = (uint64_t)60<<32 | 80, s[3].size = 251;
	s[4].info = (uint64_t)90<<32 | 100, s[4].size = 250; // not above 250
	MBT_CHECK_INT(mb_cal_high_cov(0, 5, s, 250), ==, 70); // [0,50) + [60,80)
	MBT_CHECK_INT(mb_cal_high_cov(0, 5, s, 1000), ==, 0);
	MBT_CHECK_INT(mb_cal_high_cov(0, 5, s, 400), ==, 35); // [20,50) + [70,75)
	return 0;
}

/* mb_opt_init and mb_opt_preset: the sr/adap and lr values (usage text and man page) */
static int case_presets(int argc, char *argv[])
{
	mb_opt_t o;
	mb_opt_init(&o);
	MBT_CHECK((o.flag & MB_F_PE) && (o.flag & MB_F_ADAP) && !(o.flag & MB_F_LONG));
	MBT_CHECK(o.bw == 100 && o.bw_long == 20000 && o.max_gap == 100 && o.max_sr_len == 325);
	MBT_CHECK(o.max_chain_skip == 25 && o.max_chain_iter == 5000 && o.min_chain_score == 25 && o.chain_gap_scale == 0.8f);
	MBT_CHECK(o.mask_level == 0.5f && o.mask_len == INT_MAX && o.pri_ratio == 0.5f && o.best_n == 50);
	MBT_CHECK(o.a == 2 && o.b == 8 && o.b_ts == 0 && o.b_ambi == 1 && o.q == 12 && o.q2 == 23 && o.e == 2 && o.e2 == 1);
	MBT_CHECK(o.end_bonus == 10 && o.min_dp_max == 30 && o.zdrop == 80 && o.zdrop_inv == 80 && o.min_ksw_len == 20);
	MBT_CHECK(o.max_sw_mat == 100000000 && o.min_len == 19 && o.seed == 11);
	o.a = 5, o.flag |= MB_F_EQX;
	MBT_CHECK_INT(mb_opt_preset(&o, "sr"), ==, 0); // a preset resets everything first
	MBT_CHECK_INT(o.a, ==, 2);
	MBT_CHECK((o.flag & MB_F_PE) && !(o.flag & MB_F_ADAP) && !(o.flag & MB_F_LONG) && !(o.flag & MB_F_EQX));
	MBT_CHECK(o.bw == 100 && o.max_gap == 100 && o.min_chain_score == 25 && o.min_dp_max == 30 && o.best_n == 50 && o.end_bonus == 10);
	MBT_CHECK_INT(mb_opt_preset(&o, "lr"), ==, 0);
	MBT_CHECK(!(o.flag & MB_F_PE) && !(o.flag & MB_F_ADAP) && (o.flag & MB_F_LONG));
	MBT_CHECK(o.bw == 500 && o.bw_long == 20000 && o.max_gap == 5000 && o.min_chain_score == 40 && o.best_n == 5);
	MBT_CHECK(o.end_bonus == -1 && o.min_dp_max == 50 && o.zdrop == 400 && o.zdrop_inv == 240 && o.min_ksw_len == 200);
	MBT_CHECK(o.a == 2 && o.b == 8 && o.q == 12 && o.q2 == 23 && o.e == 2 && o.e2 == 1);
	MBT_CHECK_INT(mb_opt_preset(&o, "adap"), ==, 0);
	MBT_CHECK((o.flag & MB_F_PE) && (o.flag & MB_F_ADAP) && o.bw == 100);
	MBT_CHECK_INT(mb_opt_preset(&o, "map-ont"), ==, -1);
	MBT_CHECK_INT(mb_opt_preset(&o, ""), ==, -1);
	return 0;
}

/* mb_opt_adap: x1 - (x1 - x0) * 2^(-(max(L,100)-100)/1900) for -g (5000), -w (500), -s (50),
 * -m (40), zdrop (400) and zdrop_inv (240); (x0 - 5) * 2^(...) + 5 for -N; -g at most L,
 * bw_long at most 5L but at least bw. Values computed independently for each L. */
static int case_opt_adap(int argc, char *argv[])
{
	static const int tab[][9] = { // L, max_gap, bw, bw_long, zdrop, zdrop_inv, best_n, min_dp_max, min_chain_score
		{ 50, 50, 100, 250, 80, 80, 50, 30, 25 },
		{ 100, 100, 100, 500, 80, 80, 50, 30, 25 },
		{ 150, 150, 107, 750, 86, 83, 49, 30, 25 },
		{ 250, 250, 121, 1250, 97, 89, 48, 31, 26 },
		{ 1000, 1000, 212, 5000, 170, 125, 37, 36, 29 },
		{ 2000, 2000, 300, 10000, 240, 160, 27, 40, 32 },
		{ 5000, 4180, 433, 20000, 346, 213, 13, 47, 37 },
		{ 10000, 4868, 489, 20000, 391, 236, 6, 49, 40 },
		{ 100000, 5000, 500, 20000, 400, 240, 5, 50, 40 }
	};
	mb_opt_t o0, o;
	int i;
	mb_opt_init(&o0);
	for (i = 0; i < (int)(sizeof(tab) / sizeof(tab[0])); ++i) {
		const int *t = tab[i];
		mb_opt_adap(&o0, t[0], &o);
		MBT_CHECKF(o.max_gap == t[1] && o.bw == t[2] && o.bw_long == t[3], "L=%d: max_gap %d bw %d bw_long %d", t[0], o.max_gap, o.bw, o.bw_long);
		MBT_CHECKF(o.zdrop == t[4] && o.zdrop_inv == t[5] && o.best_n == t[6], "L=%d: zdrop %d zdrop_inv %d best_n %d", t[0], o.zdrop, o.zdrop_inv, o.best_n);
		MBT_CHECKF(o.min_dp_max == t[7] && o.min_chain_score == t[8], "L=%d: min_dp_max %d min_chain_score %d", t[0], o.min_dp_max, o.min_chain_score);
		MBT_CHECKF(o.a == o0.a && o.q == o0.q && o.end_bonus == o0.end_bonus && o.min_ksw_len == o0.min_ksw_len && o.flag == o0.flag, "L=%d: other fields changed", t[0]);
	}
	// values at or above the long-read value are left alone; so is everything without MB_F_ADAP
	o0.max_gap = 6000, o0.bw = 600, o0.min_chain_score = 45, o0.best_n = 3;
	mb_opt_adap(&o0, 150, &o);
	MBT_CHECK(o.max_gap == 6000 && o.bw == 600 && o.min_chain_score == 45 && o.best_n == 3 && o.bw_long == 750);
	mb_opt_init(&o0);
	o0.flag &= ~MB_F_ADAP;
	mb_opt_adap(&o0, 5000, &o);
	MBT_CHECK(memcmp(&o, &o0, sizeof(o)) == 0);
	return 0;
}

/* mb_is_sr_mode: short-read mode unless --long, or --adap with a read longer than -l */
static int case_is_sr_mode(int argc, char *argv[])
{
	mb_opt_t o;
	mb_opt_init(&o);
	MBT_CHECK(mb_is_sr_mode(&o, 150) == 1 && mb_is_sr_mode(&o, 325) == 1 && mb_is_sr_mode(&o, 326) == 0);
	o.flag &= ~MB_F_ADAP;
	MBT_CHECK(mb_is_sr_mode(&o, 100000) == 1);
	o.flag |= MB_F_LONG;
	MBT_CHECK(mb_is_sr_mode(&o, 10) == 0);
	mb_opt_preset(&o, "lr");
	MBT_CHECK(mb_is_sr_mode(&o, 150) == 0);
	return 0;
}

/* mb_hash_str is 32-bit FNV-1a (published test vectors); mb_hash64 is the splitmix64
 * finalizer (values from an independent implementation); mb_log2 approximates log2
 * within 0.01; mb_seq_rev reverses in place. */
static int case_helpers(int argc, char *argv[])
{
	uint8_t s[5] = {0, 1, 2, 3, 4};
	double x, err = 0.0;
	MBT_CHECK(mb_hash_str("") == 0x811c9dc5U);
	MBT_CHECK(mb_hash_str("a") == 0xe40c292cU);
	MBT_CHECK(mb_hash_str("foobar") == 0xbf9cf968U);
	MBT_CHECK(mb_hash_str("read/1") == 0x82abd9bdU);
	MBT_CHECK(mb_hash64(0) == 0);
	MBT_CHECK(mb_hash64(1) == 0x5692161d100b05e5ULL);
	MBT_CHECK(mb_hash64(11) == 0x3462d848f53abb6dULL);
	MBT_CHECK(mb_hash64(0xdeadbeefULL) == 0x4e062702ec929eeaULL);
	MBT_CHECK(mb_hash64(~0ULL) == 0xb4d055fcf2cbbd7bULL);
	for (x = 2.0; x < 2e6; x *= 1.0137) {
		double d = fabs(mb_log2((float)x) - log2((float)x));
		if (d > err) err = d;
	}
	MBT_CHECKF(err < 0.01, "mb_log2 is off by %g", err);
	mb_seq_rev(5, s);
	MBT_CHECK(s[0] == 4 && s[1] == 3 && s[2] == 2 && s[3] == 1 && s[4] == 0);
	mb_seq_rev(4, s);
	MBT_CHECK(s[0] == 1 && s[3] == 4 && s[4] == 0);
	return 0;
}

static int bh_cmp64(const void *a, const void *b)
{
	uint64_t x = *(const uint64_t*)a, y = *(const uint64_t*)b;
	return x < y? -1 : x > y? 1 : 0;
}

/* radix_sort_mb64 and radix_sort_mb128x sort ascending (by .x for mb128_t) and keep
 * every element; checked against qsort */
static int case_radix_sort(int argc, char *argv[])
{
	uint64_t st = 12345, *a, *b;
	mb128_t *c;
	int n, i, k, ok;
	for (k = 0; k < 6; ++k) {
		n = k == 0? 1 : k == 1? 2 : k == 2? 63 : k == 3? 64 : k == 4? 1000 : 70000;
		a = (uint64_t*)malloc(n * 8), b = (uint64_t*)malloc(n * 8), c = (mb128_t*)malloc(n * 16);
		for (i = 0; i < n; ++i) {
			st = st * 6364136223846793005ULL + 1442695040888963407ULL;
			a[i] = k & 1? st : st >> (i % 64); // full-width keys, and keys with many leading zero bytes
			c[i].x = a[i], c[i].y = a[i] ^ 0x5555;
		}
		memcpy(b, a, n * 8);
		qsort(b, n, 8, bh_cmp64);
		radix_sort_mb64(a, a + n);
		MBT_CHECKF(memcmp(a, b, n * 8) == 0, "radix_sort_mb64 with n=%d", n);
		radix_sort_mb128x(c, c + n);
		for (i = 0, ok = 1; i < n; ++i)
			if (c[i].x != b[i] || c[i].y != (c[i].x ^ 0x5555)) ok = 0;
		MBT_CHECKF(ok, "radix_sort_mb128x with n=%d", n);
		free(a); free(b); free(c);
	}
	return 0;
}

static const mbt_case_t cases[] = {
	{ "set_parent", case_set_parent },
	{ "set_sam_pri", case_set_sam_pri },
	{ "select_sub", case_select_sub },
	{ "hit_sort", case_hit_sort },
	{ "filter_hits", case_filter_hits },
	{ "set_mapq", case_set_mapq },
	{ "high_cov", case_high_cov },
	{ "presets", case_presets },
	{ "opt_adap", case_opt_adap },
	{ "is_sr_mode", case_is_sr_mode },
	{ "helpers", case_helpers },
	{ "radix_sort", case_radix_sort },
	{ 0, 0 }
};

MBT_MAIN(cases)
