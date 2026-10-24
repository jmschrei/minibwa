/* Area B driver: the ksw2 kernels (ksw_extz2_sse, ksw_extd2_sse, ksw_ll_i16,
 * ksw_ll_u8_core) and the inline helpers in ksw2.h, checked against a
 * brute-force affine-gap dynamic program written here.
 *
 * Conventions (ksw2.h): sequences are 0..4 codes (4 = N), mat is a 5x5
 * matrix indexed [target * 5 + query], a gap of length l costs q + l*e and,
 * in the dual-affine kernel, min(q + l*e, q2 + l*e2). Every random input
 * comes from a splitmix64 stream seeded from the command line. */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "ksw2.h"
#include "minibwa.h"
#include "minibwa_test.h"

#define BF_NEG (-0x20000000)

static uint64_t bk_state;

static uint64_t bk_rand(void)
{
	uint64_t z = (bk_state += 0x9E3779B97F4A7C15ULL);
	z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
	z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
	return z ^ (z >> 31);
}

static int bk_below(int n) { return (int)(bk_rand() % (uint64_t)n); }

static inline int bk_max2(int a, int b) { return a > b? a : b; }

/* the standard nt4 matrix built by hand: a on the diagonal, -b off it, -ambi in the N row and column */
static void bk_mat(int8_t mat[25], int a, int b, int b_ts, int ambi)
{
	int i, j;
	for (i = 0; i < 5; ++i)
		for (j = 0; j < 5; ++j)
			mat[i * 5 + j] = i == 4 || j == 4? -ambi : i == j? a : (b_ts && (i ^ j) == 2)? -b_ts : -b;
}

/* a random pair: target of length tlen, query derived from it (or unrelated) with
 * substitutions, indels and N at the given per-mille rates */
static void bk_pair(int tlen, uint8_t *t, int *qlen, uint8_t *q, int qmax, int sub_pm, int indel_pm, int n_pm, int unrelated)
{
	int i, k = 0;
	for (i = 0; i < tlen; ++i)
		t[i] = bk_below(1000) < n_pm? 4 : bk_below(4);
	if (unrelated) {
		*qlen = 1 + bk_below(qmax);
		for (i = 0; i < *qlen; ++i) q[i] = bk_below(1000) < n_pm? 4 : bk_below(4);
		return;
	}
	for (i = 0; i < tlen && k < qmax; ++i) {
		int x = bk_below(1000);
		if (x < indel_pm / 2) { // deletion of 1-6 bases
			i += bk_below(6);
			continue;
		} else if (x < indel_pm) { // insertion of 1-6 bases
			int l = 1 + bk_below(6), j;
			for (j = 0; j < l && k < qmax; ++j) q[k++] = bk_below(4);
		}
		if (k >= qmax) break;
		x = bk_below(1000);
		q[k++] = x < n_pm? 4 : x < n_pm + sub_pm? (t[i] + 1 + bk_below(3)) % 4 : t[i];
	}
	if (k == 0) q[k++] = bk_below(4);
	*qlen = k;
}

/* Brute-force Gotoh with two gap models. H, E1, E2, F1, F2 are (tlen+1)*(qlen+1);
 * row i = target prefix length, column j = query prefix length. E: gap in the query
 * (CIGAR D, consumes target); F: gap in the target (CIGAR I, consumes query).
 * local: scores floored at 0 (Smith-Waterman). */
typedef struct {
	int tlen, qlen;
	int *H;
} bk_dp_t;

static inline int bk_H(const bk_dp_t *d, int i, int j) { return d->H[i * (d->qlen + 1) + j]; }

static void bk_dp(bk_dp_t *d, int qlen, const uint8_t *q, int tlen, const uint8_t *t, const int8_t *mat, int q1, int e1, int q2, int e2, int local)
{
	int i, j, n = (tlen + 1) * (qlen + 1);
	int *E1 = (int*)malloc(n * sizeof(int)), *E2 = (int*)malloc(n * sizeof(int));
	int *F1 = (int*)malloc(n * sizeof(int)), *F2 = (int*)malloc(n * sizeof(int));
	d->tlen = tlen, d->qlen = qlen;
	d->H = (int*)malloc(n * sizeof(int));
#define IDX(i, j) ((i) * (qlen + 1) + (j))
	for (i = 0; i <= tlen; ++i) {
		for (j = 0; j <= qlen; ++j) {
			int h, e1v = BF_NEG, e2v = BF_NEG, f1v = BF_NEG, f2v = BF_NEG;
			if (i > 0) {
				e1v = bk_max2(d->H[IDX(i-1, j)] - (q1 + e1), E1[IDX(i-1, j)] - e1);
				e2v = bk_max2(d->H[IDX(i-1, j)] - (q2 + e2), E2[IDX(i-1, j)] - e2);
			}
			if (j > 0) {
				f1v = bk_max2(d->H[IDX(i, j-1)] - (q1 + e1), F1[IDX(i, j-1)] - e1);
				f2v = bk_max2(d->H[IDX(i, j-1)] - (q2 + e2), F2[IDX(i, j-1)] - e2);
			}
			if (i == 0 && j == 0) h = 0;
			else {
				h = BF_NEG;
				if (i > 0 && j > 0) h = d->H[IDX(i-1, j-1)] + mat[t[i-1] * 5 + q[j-1]];
				h = bk_max2(h, bk_max2(bk_max2(e1v, e2v), bk_max2(f1v, f2v)));
				if (local && h < 0) h = 0;
			}
			if (local && i == 0) h = 0;
			if (local && j == 0) h = 0;
			d->H[IDX(i, j)] = h;
			E1[IDX(i, j)] = e1v, E2[IDX(i, j)] = e2v, F1[IDX(i, j)] = f1v, F2[IDX(i, j)] = f2v;
		}
	}
#undef IDX
	free(E1); free(E2); free(F1); free(F2);
}

/* score of a CIGAR from (0,0); spans returned in *qs and *ts; -1 on an unknown op */
static int bk_cigar_score(int n_cigar, const uint32_t *cigar, const uint8_t *q, const uint8_t *t, const int8_t *mat, int q1, int e1, int q2, int e2, int *qs, int *ts)
{
	int k, i = 0, j = 0, s = 0;
	for (k = 0; k < n_cigar; ++k) {
		int op = cigar[k] & 0xf, len = cigar[k] >> 4, l;
		if (op == KSW_CIGAR_MATCH) {
			for (l = 0; l < len; ++l) s += mat[t[i + l] * 5 + q[j + l]];
			i += len, j += len;
		} else if (op == KSW_CIGAR_INS) {
			s -= q1 + len * e1 < q2 + len * e2? q1 + len * e1 : q2 + len * e2;
			j += len;
		} else if (op == KSW_CIGAR_DEL) {
			s -= q1 + len * e1 < q2 + len * e2? q1 + len * e1 : q2 + len * e2;
			i += len;
		} else return -1;
	}
	*qs = j, *ts = i;
	return s;
}

static int bk_cigar_ok(int n_cigar, const uint32_t *cigar)
{
	int k;
	for (k = 0; k < n_cigar; ++k) {
		if ((cigar[k] >> 4) == 0) return 0; // no zero-length operation
		if (k > 0 && (cigar[k] & 0xf) == (cigar[k-1] & 0xf)) return 0; // adjacent operations differ
	}
	return 1;
}

static void bk_rev_cigar(int n, uint32_t *c)
{
	int i;
	for (i = 0; i < n>>1; ++i) {
		uint32_t x = c[i];
		c[i] = c[n-1-i], c[n-1-i] = x;
	}
}

typedef struct {
	int dual;    // ksw_extd2_sse instead of ksw_extz2_sse
	int q, e, q2, e2;
} bk_gap_t;

static void bk_run(const bk_gap_t *g, int qlen, const uint8_t *q, int tlen, const uint8_t *t, const int8_t *mat, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez)
{
	if (g->dual) ksw_extd2_sse(0, qlen, q, tlen, t, 5, mat, g->q, g->e, g->q2, g->e2, w, zdrop, end_bonus, flag, ez);
	else ksw_extz2_sse(0, qlen, q, tlen, t, 5, mat, g->q, g->e, w, zdrop, end_bonus, flag, ez);
}

/* Global alignment (no KSW_EZ_EXTZ_ONLY), full band, no Z-drop: ez->score is the
 * optimal global score and the CIGAR spans both sequences with that score.
 * args: seed n_iter dual [generic [q e q2 e2]] (default gap costs 12,2 and 23,1) */
static int case_global(int argc, char *argv[])
{
	bk_gap_t g;
	int it, n_iter, generic;
	uint8_t t[200], q[200];
	int8_t mat[25];
	ksw_extz_t ez;
	if (argc < 3) return MBT_FAIL("usage: global <seed> <n_iter> <dual> [generic [q e q2 e2]]"), 1;
	bk_state = strtoull(argv[0], 0, 10);
	n_iter = atoi(argv[1]);
	g.dual = atoi(argv[2]);
	generic = argc > 3? atoi(argv[3]) : 0;
	g.q = 12, g.e = 2, g.q2 = g.dual? 23 : 12, g.e2 = g.dual? 1 : 2;
	if (argc > 7) g.q = atoi(argv[4]), g.e = atoi(argv[5]), g.q2 = atoi(argv[6]), g.e2 = atoi(argv[7]);
	bk_mat(mat, 2, 8, generic? 3 : 0, 1);
	memset(&ez, 0, sizeof(ez));
	for (it = 0; it < n_iter; ++it) {
		int tlen = 1 + bk_below(90), qlen, flag, qs, ts, s, bf;
		bk_dp_t d;
		bk_pair(tlen, t, &qlen, q, 120, 60, 40, 10, bk_below(5) == 0);
		bk_dp(&d, qlen, q, tlen, t, mat, g.q, g.e, g.q2, g.e2, 0);
		bf = bk_H(&d, tlen, qlen);
		flag = (it & 1? KSW_EZ_RIGHT : 0) | (generic? KSW_EZ_GENERIC_SC : 0);
		bk_run(&g, qlen, q, tlen, t, mat, -1, -1, 0, flag, &ez);
		MBT_CHECKF(ez.score == bf, "iter %d (tlen %d qlen %d flag %x): score %d, brute force %d", it, tlen, qlen, flag, ez.score, bf);
		MBT_CHECKF(!ez.zdropped, "iter %d: Z-dropped with zdrop disabled", it);
		s = bk_cigar_score(ez.n_cigar, ez.cigar, q, t, mat, g.q, g.e, g.q2, g.e2, &qs, &ts);
		MBT_CHECKF(qs == qlen && ts == tlen, "iter %d: CIGAR spans %d/%d, sequences %d/%d", it, qs, ts, qlen, tlen);
		MBT_CHECKF(s == bf, "iter %d: CIGAR scores %d, optimum %d", it, s, bf);
		MBT_CHECKF(bk_cigar_ok(ez.n_cigar, ez.cigar), "iter %d: CIGAR has zero-length or unmerged operations", it);
		// the score-only mode agrees
		bk_run(&g, qlen, q, tlen, t, mat, -1, -1, 0, flag | KSW_EZ_SCORE_ONLY, &ez);
		MBT_CHECKF(ez.score == bf && ez.n_cigar == 0, "iter %d: score-only mode gives %d (n_cigar %d), optimum %d", it, ez.score, ez.n_cigar, bf);
		free(d.H);
	}
	free(ez.cigar);
	return 0;
}

/* Extension (KSW_EZ_EXTZ_ONLY), full band, no Z-drop: max is the best score of any
 * prefix pair (at least 0); the CIGAR goes to the query end, scoring the best
 * end-reaching score mqe, iff mqe + end_bonus > max, and to the max cell otherwise.
 * args: seed n_iter dual */
static int case_extension(int argc, char *argv[])
{
	bk_gap_t g;
	int it, n_iter;
	uint8_t t[200], q[200];
	int8_t mat[25];
	ksw_extz_t ez;
	uint32_t *c2 = 0;
	if (argc < 3) return MBT_FAIL("usage: extension <seed> <n_iter> <dual>"), 1;
	bk_state = strtoull(argv[0], 0, 10);
	n_iter = atoi(argv[1]);
	g.dual = atoi(argv[2]);
	g.q = 12, g.e = 2, g.q2 = g.dual? 23 : 12, g.e2 = g.dual? 1 : 2;
	bk_mat(mat, 2, 8, 0, 1);
	memset(&ez, 0, sizeof(ez));
	for (it = 0; it < n_iter; ++it) {
		int tlen = 1 + bk_below(90), qlen, i, j, bmax = 0, bmqe = KSW_NEG_INF, bmte = KSW_NEG_INF, end_bonus, qs, ts, s, n2;
		bk_dp_t d;
		bk_pair(tlen, t, &qlen, q, 120, 60, 40, 10, bk_below(4) == 0);
		bk_dp(&d, qlen, q, tlen, t, mat, g.q, g.e, g.q2, g.e2, 0);
		for (i = 1; i <= tlen; ++i)
			for (j = 1; j <= qlen; ++j)
				bmax = bk_max2(bmax, bk_H(&d, i, j));
		for (i = 1; i <= tlen; ++i) bmqe = bk_max2(bmqe, bk_H(&d, i, qlen));
		for (j = 1; j <= qlen; ++j) bmte = bk_max2(bmte, bk_H(&d, tlen, j));
		end_bonus = bk_below(3) * 5;
		bk_run(&g, qlen, q, tlen, t, mat, -1, -1, end_bonus, KSW_EZ_EXTZ_ONLY | (it & 1? KSW_EZ_RIGHT : 0), &ez);
		MBT_CHECKF((int)ez.max == bmax, "iter %d: max %d, brute force %d", it, ez.max, bmax);
		if (bmax > 0)
			MBT_CHECKF(ez.max_t >= 0 && ez.max_q >= 0 && bk_H(&d, ez.max_t + 1, ez.max_q + 1) == bmax, "iter %d: max cell (%d,%d) does not score %d", it, ez.max_t, ez.max_q, bmax);
		// mqe/mte: the best end-reaching scores seen; exact whenever the CIGAR goes to the
		// end (ksw_extd2_sse stops early once no cell can beat max, so not otherwise)
		MBT_CHECKF(ez.reach_end == (bmqe + end_bonus > bmax), "iter %d: reach_end %d with mqe %d + bonus %d vs max %d", it, ez.reach_end, bmqe, end_bonus, bmax);
		MBT_CHECKF(ez.mqe <= bmqe && (ez.mqe_t < 0 || bk_H(&d, ez.mqe_t + 1, qlen) == ez.mqe), "iter %d: mqe %d at t=%d, brute force %d", it, ez.mqe, ez.mqe_t, bmqe);
		MBT_CHECKF(ez.mte <= bmte && (ez.mte_q < 0 || bk_H(&d, tlen, ez.mte_q + 1) == ez.mte), "iter %d: mte %d at q=%d, brute force %d", it, ez.mte, ez.mte_q, bmte);
		if (ez.reach_end) MBT_CHECKF(ez.mqe == bmqe, "iter %d: reaches the end with mqe %d, brute force %d", it, ez.mqe, bmqe);
		s = bk_cigar_score(ez.n_cigar, ez.cigar, q, t, mat, g.q, g.e, g.q2, g.e2, &qs, &ts);
		if (ez.reach_end) {
			MBT_CHECKF(qs == qlen && ts == ez.mqe_t + 1 && s == bmqe, "iter %d: end-reaching CIGAR spans %d/%d scoring %d; want %d/%d scoring %d", it, qs, ts, s, qlen, ez.mqe_t + 1, bmqe);
		} else if (bmax > 0) {
			MBT_CHECKF(qs == ez.max_q + 1 && ts == ez.max_t + 1 && s == bmax, "iter %d: CIGAR spans %d/%d scoring %d; want %d/%d scoring %d", it, qs, ts, s, ez.max_q + 1, ez.max_t + 1, bmax);
		}
		MBT_CHECKF(bk_cigar_ok(ez.n_cigar, ez.cigar), "iter %d: CIGAR has zero-length or unmerged operations", it);
		// KSW_EZ_REV_CIGAR returns the same CIGAR reversed
		n2 = ez.n_cigar;
		c2 = (uint32_t*)realloc(c2, (n2 + 1) * 4);
		memcpy(c2, ez.cigar, n2 * 4);
		bk_run(&g, qlen, q, tlen, t, mat, -1, -1, end_bonus, KSW_EZ_EXTZ_ONLY | KSW_EZ_REV_CIGAR | (it & 1? KSW_EZ_RIGHT : 0), &ez);
		bk_rev_cigar(ez.n_cigar, ez.cigar);
		MBT_CHECKF(ez.n_cigar == n2 && memcmp(ez.cigar, c2, n2 * 4) == 0, "iter %d: KSW_EZ_REV_CIGAR is not the reversed CIGAR", it);
		free(d.H);
	}
	free(c2);
	free(ez.cigar);
	return 0;
}

/* KSW_EZ_APPROX_MAX in global mode (how mb_align1 fills gaps between anchors): the
 * global score is still exact and agrees with the CIGAR. args: seed n_iter dual */
static int case_global_approx(int argc, char *argv[])
{
	bk_gap_t g;
	int it, n_iter;
	uint8_t t[200], q[200];
	int8_t mat[25];
	ksw_extz_t ez;
	if (argc < 3) return MBT_FAIL("usage: global_approx <seed> <n_iter> <dual>"), 1;
	bk_state = strtoull(argv[0], 0, 10);
	n_iter = atoi(argv[1]);
	g.dual = atoi(argv[2]);
	g.q = 12, g.e = 2, g.q2 = g.dual? 23 : 12, g.e2 = g.dual? 1 : 2;
	bk_mat(mat, 2, 8, 0, 1);
	memset(&ez, 0, sizeof(ez));
	for (it = 0; it < n_iter; ++it) {
		int tlen = 1 + bk_below(90), qlen, qs, ts, s, bf;
		bk_dp_t d;
		bk_pair(tlen, t, &qlen, q, 120, 40, 30, 5, 0);
		bk_dp(&d, qlen, q, tlen, t, mat, g.q, g.e, g.q2, g.e2, 0);
		bf = bk_H(&d, tlen, qlen);
		bk_run(&g, qlen, q, tlen, t, mat, -1, -1, -1, KSW_EZ_APPROX_MAX, &ez);
		s = bk_cigar_score(ez.n_cigar, ez.cigar, q, t, mat, g.q, g.e, g.q2, g.e2, &qs, &ts);
		MBT_CHECKF(ez.score == bf, "iter %d: approx-max global score %d, optimum %d", it, ez.score, bf);
		MBT_CHECKF(qs == qlen && ts == tlen && s == bf, "iter %d: CIGAR spans %d/%d scoring %d", it, qs, ts, s);
		free(d.H);
	}
	free(ez.cigar);
	return 0;
}

/* A band of half-width w (ksw2 computes cells with 1-w <= t-q <= w): when the
 * unbanded optimal path stays at least `margin` cells inside the band, the banded run
 * finds the optimum and its CIGAR scores it. Nearer the band edge ksw2 is an
 * approximation (its out-of-band neighbours are not -inf), so nothing is asserted
 * there. args: seed n_iter dual margin */
static int case_band(int argc, char *argv[])
{
	bk_gap_t g;
	int it, n_iter, n_fit = 0, margin;
	uint8_t t[200], q[200];
	int8_t mat[25];
	ksw_extz_t ez;
	if (argc < 4) return MBT_FAIL("usage: band <seed> <n_iter> <dual> <margin>"), 1;
	bk_state = strtoull(argv[0], 0, 10);
	n_iter = atoi(argv[1]);
	g.dual = atoi(argv[2]);
	margin = atoi(argv[3]);
	g.q = 12, g.e = 2, g.q2 = g.dual? 23 : 12, g.e2 = g.dual? 1 : 2;
	bk_mat(mat, 2, 8, 0, 1);
	memset(&ez, 0, sizeof(ez));
	for (it = 0; it < n_iter; ++it) {
		int tlen = 20 + bk_below(80), qlen, qs, ts, s, bf, w, k, i, j, max_off = 0, min_off = 0;
		bk_dp_t d;
		bk_pair(tlen, t, &qlen, q, 120, 30, 30, 0, 0);
		bk_dp(&d, qlen, q, tlen, t, mat, g.q, g.e, g.q2, g.e2, 0);
		bf = bk_H(&d, tlen, qlen);
		free(d.H);
		bk_run(&g, qlen, q, tlen, t, mat, -1, -1, 0, 0, &ez); // the full-band path
		for (k = i = j = 0; k < ez.n_cigar; ++k) {
			int op = ez.cigar[k] & 0xf, len = ez.cigar[k] >> 4;
			if (op == KSW_CIGAR_MATCH) i += len, j += len;
			else if (op == KSW_CIGAR_DEL) i += len;
			else j += len;
			if (i - j > max_off) max_off = i - j;
			if (i - j < min_off) min_off = i - j;
		}
		w = (max_off > -min_off? max_off : -min_off) + margin + bk_below(4);
		if (w >= (tlen > qlen? tlen : qlen)) continue; // not a band
		++n_fit;
		bk_run(&g, qlen, q, tlen, t, mat, w, -1, 0, 0, &ez);
		s = bk_cigar_score(ez.n_cigar, ez.cigar, q, t, mat, g.q, g.e, g.q2, g.e2, &qs, &ts);
		MBT_CHECKF(!ez.zdropped && ez.score == bf, "iter %d: path offsets %d..%d, band %d: score %d, optimum %d", it, min_off, max_off, w, ez.score, bf);
		MBT_CHECKF(qs == qlen && ts == tlen && s == bf, "iter %d: band %d CIGAR spans %d/%d scoring %d, optimum %d", it, w, qs, ts, s, bf);
	}
	MBT_CHECKF(n_fit >= n_iter / 2, "only %d of %d runs used a band narrower than the sequences", n_fit, n_iter);
	free(ez.cigar);
	return 0;
}

/* Z-drop: 100 identical bases followed by 100 bases that mismatch on the diagonal.
 * The best extension is the 100 matches (score 200 at (99,99)); a Z-drop threshold of
 * 50 stops the extension, -1 and a huge threshold do not. args: seed dual */
static int case_zdrop(int argc, char *argv[])
{
	bk_gap_t g;
	uint8_t t[200], q[200];
	int8_t mat[25];
	int i, z;
	ksw_extz_t ez;
	if (argc < 2) return MBT_FAIL("usage: zdrop <seed> <dual>"), 1;
	bk_state = strtoull(argv[0], 0, 10);
	g.dual = atoi(argv[1]);
	g.q = 12, g.e = 2, g.q2 = g.dual? 23 : 12, g.e2 = g.dual? 1 : 2;
	bk_mat(mat, 2, 8, 0, 1);
	memset(&ez, 0, sizeof(ez));
	for (i = 0; i < 100; ++i) t[i] = q[i] = bk_below(4);
	for (i = 100; i < 200; ++i) t[i] = bk_below(4), q[i] = (t[i] + 1 + bk_below(3)) % 4;
	{
		bk_dp_t d;
		int j, bmax = 0;
		bk_dp(&d, 200, q, 200, t, mat, g.q, g.e, g.q2, g.e2, 0);
		for (i = 1; i <= 200; ++i)
			for (j = 1; j <= 200; ++j)
				bmax = bk_max2(bmax, bk_H(&d, i, j));
		MBT_REQUIRE(bmax == 200 && bk_H(&d, 100, 100) == 200);
		free(d.H);
	}
	for (z = 0; z < 3; ++z) {
		int zdrop = z == 0? 50 : z == 1? -1 : 100000;
		bk_run(&g, 200, q, 200, t, mat, -1, zdrop, 0, KSW_EZ_EXTZ_ONLY, &ez);
		MBT_CHECKF(ez.zdropped == (zdrop == 50), "zdrop %d: zdropped=%d", zdrop, ez.zdropped);
		MBT_CHECKF(ez.max == 200 && ez.max_t == 99 && ez.max_q == 99, "zdrop %d: max %d at (%d,%d)", zdrop, ez.max, ez.max_t, ez.max_q);
		MBT_CHECKF(ez.n_cigar == 1 && ez.cigar[0] == (100<<4 | KSW_CIGAR_MATCH), "zdrop %d: CIGAR is not 100M", zdrop);
		// global mode: a Z-drop leaves no global score; the CIGAR runs to the best cell
		bk_run(&g, 200, q, 200, t, mat, -1, zdrop, 0, 0, &ez);
		if (zdrop == 50) MBT_CHECKF(ez.zdropped && ez.score == KSW_NEG_INF && ez.max == 200 && ez.n_cigar == 1 && ez.cigar[0] == (100<<4 | KSW_CIGAR_MATCH),
				"global zdrop: zdropped=%d score=%d max=%d n_cigar=%d", ez.zdropped, ez.score, ez.max, ez.n_cigar);
		else MBT_CHECKF(!ez.zdropped && ez.score > KSW_NEG_INF, "global without zdrop: zdropped=%d score=%d", ez.zdropped, ez.score);
	}
	free(ez.cigar);
	return 0;
}

/* Local alignment (ksw_ll_i16, and ksw_ll_u8_core with 8-bit cells) against
 * brute-force Smith-Waterman with gap cost q + l*e. args: seed n_iter */
static int case_local(int argc, char *argv[])
{
	int it, n_iter;
	uint8_t t[300], q[200];
	int8_t mat[25];
	if (argc < 2) return MBT_FAIL("usage: local <seed> <n_iter>"), 1;
	bk_state = strtoull(argv[0], 0, 10);
	n_iter = atoi(argv[1]);
	bk_mat(mat, 2, 8, 0, 1);
	for (it = 0; it < n_iter; ++it) {
		int tlen, qlen, i, j, bf = 0, sc, qe, te, k, flank = bk_below(60);
		bk_dp_t d;
		void *qp;
		ksw_llrst_t r;
		// a query that matches somewhere inside a longer target
		tlen = 1 + bk_below(80);
		bk_pair(tlen, t + flank, &qlen, q, 100, 80, 40, 10, bk_below(5) == 0);
		for (i = 0; i < flank; ++i) t[i] = bk_below(4);
		k = bk_below(60);
		for (i = 0; i < k; ++i) t[flank + tlen + i] = bk_below(4);
		tlen += flank + k;
		bk_dp(&d, qlen, q, tlen, t, mat, 12, 2, 12, 2, 1);
		for (i = 1; i <= tlen; ++i)
			for (j = 1; j <= qlen; ++j)
				bf = bk_max2(bf, bk_H(&d, i, j));
		qp = ksw_ll_qinit(0, 2, qlen, q, 5, mat);
		sc = ksw_ll_i16(qp, tlen, t, 12, 2, &qe, &te);
		free(qp);
		MBT_CHECKF(sc == bf, "iter %d: ksw_ll_i16 score %d, brute force %d", it, sc, bf);
		if (bf > 0) MBT_CHECKF(te >= 0 && qe >= 0 && te < tlen && qe < qlen && bk_H(&d, te + 1, qe + 1) == bf, "iter %d: end (%d,%d) does not score %d", it, te, qe, bf);
		qp = ksw_ll_qinit(0, 1, qlen, q, 5, mat);
		r = ksw_ll_u8_core(qp, tlen, t, 12, 2, 0);
		free(qp);
		if (bf + 8 < 255) {
			MBT_CHECKF(r.score == bf, "iter %d: ksw_ll_u8_core score %d, brute force %d", it, r.score, bf);
			if (bf > 0) MBT_CHECKF(r.te >= 0 && r.qe >= 0 && bk_H(&d, r.te + 1, r.qe + 1) == bf, "iter %d: u8 end (%d,%d) does not score %d", it, r.te, r.qe, bf);
		} else MBT_CHECKF(r.score == 255, "iter %d: u8 score %d did not saturate (brute force %d)", it, r.score, bf);
		free(d.H);
	}
	return 0;
}

/* ksw_gen_nt4_mat: a on the diagonal, -b off it, -b_ambi for N, -b_ts for transitions,
 * the bisulfite columns for mt 1 (reference C, read T) and 2 (reference G, read A);
 * signs of the inputs are normalized. */
static int case_gen_mat(int argc, char *argv[])
{
	int8_t mat[25], want[25];
	int i;
	ksw_gen_nt4_mat(mat, 2, 8, 0, 1, 0);
	bk_mat(want, 2, 8, 0, 1);
	MBT_CHECK(memcmp(mat, want, 25) == 0);
	ksw_gen_nt4_mat(mat, -2, -8, 0, -1, 0); // signs do not matter
	MBT_CHECK(memcmp(mat, want, 25) == 0);
	bk_mat(want, 1, 4, 0, 1);
	ksw_gen_nt4_mat(mat, 1, 4, 4, 1, 0); // b_ts equal to b changes nothing
	MBT_CHECK(memcmp(mat, want, 25) == 0);
	ksw_gen_nt4_mat(mat, 2, 8, 3, 2, 0); // transitions A<->G, C<->T
	bk_mat(want, 2, 8, 3, 2);
	MBT_CHECK(memcmp(mat, want, 25) == 0);
	MBT_CHECK_INT(mat[0 * 5 + 2], ==, -3);
	MBT_CHECK_INT(mat[3 * 5 + 1], ==, -3);
	MBT_CHECK_INT(mat[0 * 5 + 1], ==, -8);
	MBT_CHECK_INT(mat[4 * 5 + 4], ==, -2);
	ksw_gen_nt4_mat(mat, 2, 8, 0, 1, 1); // C2T: reference C read as T scores a match
	bk_mat(want, 2, 8, 0, 1);
	want[1 * 5 + 3] = 2;
	MBT_CHECK(memcmp(mat, want, 25) == 0);
	MBT_CHECK_INT(mat[3 * 5 + 1], ==, -8); // not symmetric
	ksw_gen_nt4_mat(mat, 2, 8, 0, 1, 2); // G2A
	bk_mat(want, 2, 8, 0, 1);
	want[2 * 5 + 0] = 2;
	MBT_CHECK(memcmp(mat, want, 25) == 0);
	for (i = 0; i < 25; ++i)
		if (mat[i] != want[i]) MBT_FAIL("G2A cell %d: %d vs %d", i, mat[i], want[i]);
	return 0;
}

/* ksw_push_cigar merges a run with the previous operation of the same type and grows
 * the array as needed. */
static int case_push_cigar(int argc, char *argv[])
{
	uint32_t *c = 0;
	int n = 0, m = 0, i;
	c = ksw_push_cigar(0, &n, &m, c, KSW_CIGAR_MATCH, 5);
	c = ksw_push_cigar(0, &n, &m, c, KSW_CIGAR_MATCH, 3);
	c = ksw_push_cigar(0, &n, &m, c, KSW_CIGAR_INS, 2);
	c = ksw_push_cigar(0, &n, &m, c, KSW_CIGAR_DEL, 1);
	c = ksw_push_cigar(0, &n, &m, c, KSW_CIGAR_DEL, 4);
	MBT_REQUIRE(n == 3);
	MBT_CHECK(c[0] == (8<<4 | KSW_CIGAR_MATCH) && c[1] == (2<<4 | KSW_CIGAR_INS) && c[2] == (5<<4 | KSW_CIGAR_DEL));
	for (i = 0; i < 1000; ++i)
		c = ksw_push_cigar(0, &n, &m, c, i & 1? KSW_CIGAR_MATCH : KSW_CIGAR_INS, 1 + i % 7);
	MBT_CHECK_INT(n, ==, 1003);
	MBT_CHECK(m >= n);
	MBT_CHECK(c[1002] == ((1 + 999 % 7) << 4 | KSW_CIGAR_MATCH));
	free(c);
	return 0;
}

/* a scoring matrix whose worst score exceeds twice q+e makes ksw_extz2_sse return no
 * alignment (its 8-bit differences cannot hold it) */
static int case_extz2_limits(int argc, char *argv[])
{
	uint8_t t[10] = {0,1,2,3,0,1,2,3,0,1}, q[10] = {0,1,2,3,0,1,2,3,0,1};
	int8_t mat[25];
	ksw_extz_t ez;
	memset(&ez, 0, sizeof(ez));
	bk_mat(mat, 2, 8, 0, 1);
	ksw_extz2_sse(0, 10, q, 10, t, 5, mat, 12, 2, -1, -1, 0, 0, &ez);
	MBT_CHECK(ez.score == 20 && ez.n_cigar == 1);
	bk_mat(mat, 2, 30, 0, 1);
	ksw_extz2_sse(0, 10, q, 10, t, 5, mat, 12, 2, -1, -1, 0, 0, &ez);
	MBT_CHECK(ez.score == KSW_NEG_INF && ez.n_cigar == 0);
	ksw_extz2_sse(0, 0, q, 10, t, 5, mat, 12, 2, -1, -1, 0, 0, &ez); // empty query
	MBT_CHECK(ez.score == KSW_NEG_INF && ez.n_cigar == 0 && ez.max == 0);
	free(ez.cigar);
	return 0;
}

static const mbt_case_t cases[] = {
	{ "global", case_global },
	{ "extension", case_extension },
	{ "global_approx", case_global_approx },
	{ "band", case_band },
	{ "zdrop", case_zdrop },
	{ "local", case_local },
	{ "gen_mat", case_gen_mat },
	{ "push_cigar", case_push_cigar },
	{ "extz2_limits", case_extz2_limits },
	{ 0, 0 }
};

MBT_MAIN(cases)
