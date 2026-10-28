#ifndef KSW_EXTD2_BODY
#include <string.h>
#include <assert.h>
#include "ksw2.h"

#if defined(__ARM_NEON)
#define __SSE4_1__
#include "s2n-lite.h"
#elif defined(__SSE4_1__)
#include <smmintrin.h>
#elif defined(__SSE2__)
#include <xmmintrin.h>
#else
#error "Missing SSE2 or NEON intrinsics"
#endif

static inline __m128i ksw_i8x4_to_i32x4(const int8_t *x)
{
#if defined(__ARM_NEON)
    return vreinterpretq_u8_s32(vmovl_s16(vget_low_s16(vmovl_s8(vcreate_s8((uint64_t)(uint32_t)*(int32_t*)x)))));
#elif defined(__SSE4_1__)
	return _mm_cvtepi8_epi32(_mm_cvtsi32_si128(*(int32_t*)x));
#else
	return _mm_setr_epi32(x[0], x[1], x[2], x[3]);
#endif
}

/* The exact-max block keeps H[] in int16 (8 cells per SSE op) on x86; s2n-lite.h has no 16-bit add
 * or widening, so NEON keeps the int32 array. The int16 array is exact while every H value written
 * lies within +-KSW_H16_LIM: a new value is an old one plus an int8 difference, so it still fits in
 * int16. After each anti-diagonal whose values leave that range, H moves to the int32 array for the
 * rest of the call. KSW_NO_H16 and a smaller KSW_H16_LIM force the int32 paths, for testing. */
#if !defined(__ARM_NEON) && !defined(KSW_NO_H16)
#define KSW_H16 1
#ifndef KSW_H16_LIM
#define KSW_H16_LIM 32640  /* 32767 - 127 and -32768 + 128 */
#endif
#define KSW_H16_MAX_T 0x7ff0 /* positions are int16 too: tlen below this */

static inline __m128i ksw_i8x8_to_i16x8(__m128i x) /* sign-extend the low 8 bytes */
{
#if defined(__SSE4_1__)
	return _mm_cvtepi8_epi16(x);
#else
	return _mm_srai_epi16(_mm_unpacklo_epi8(x, x), 8);
#endif
}

static inline __m128i ksw_blend16(__m128i a, __m128i b, __m128i mask) /* mask? b : a */
{
#if defined(__SSE4_1__)
	return _mm_blendv_epi8(a, b, mask);
#else
	return _mm_or_si128(_mm_and_si128(mask, b), _mm_andnot_si128(mask, a));
#endif
}
#endif

/* Whether ksw_apply_zdrop(ez, 1, max_H, r, max_t, zdrop, e) reads max_t: it stores it when max_H > ez->max, and otherwise
 * uses it only in a Z-drop test that cannot pass unless ez->max - max_H > zdrop (zd_all covers e < 0 and overflow). When
 * this is 0, that call returns 0 and changes nothing, so ksw_extd2_sse finds max_t and makes the call only when it is 1. */
static inline int ksw_extd2_need_max_t(const ksw_extz_t *ez, int32_t max_H, int zdrop, int zd_all)
{
	return zd_all || max_H > (int32_t)ez->max || (zdrop >= 0 && (int32_t)ez->max - max_H > zdrop);
}

#ifndef KSW_BOUND_GAP
#define KSW_BOUND_GAP 64 // KSW_EZ_EXTZ_BOUND: the most anti-diagonals between two tries of the bound
#endif
#ifndef KSW_BOUND_DROP
#define KSW_BOUND_DROP 2 // and the fall per anti-diagonal it assumes when spacing them; the result does not depend on either
#endif

/* KSW_EZ_EXTZ_BOUND: the largest H + m * (matches left) over the in-band cells of anti-diagonals r and r-1. A cell (r,t) has
 * a + t query and b - t target bases left, with a = qlen-1-r and b = tlen-1, and every cell's H is at most its diagonal
 * predecessor's plus m, so this bounds every cell those cells lead to by diagonal moves. Cells of r are t in [st0, en0], at
 * H[t]; cells of r-1 are t in [lo1, hi1] (lo1 >= st0-1, hi1 <= en0), at H[t] - v8[t] (the update this anti-diagonal made),
 * except t = st0-1, which this anti-diagonal did not update. In 32-bit lanes, exact while H + m * (a + b + 1) fits. */
#define KSW_EXTD2_MIN32(x, y, c) (c = _mm_cmpgt_epi32(x, y), _mm_or_si128(_mm_and_si128(c, y), _mm_andnot_si128(c, x)))
#define KSW_EXTD2_MAX32(x, y, c) (c = _mm_cmpgt_epi32(x, y), _mm_or_si128(_mm_and_si128(c, x), _mm_andnot_si128(c, y)))
#define KSW_EXTD2_BOUND_BODY(H) \
	int32_t t = st0, M = KSW_NEG_INF, en = hi1 < en0? hi1 : en0, end = st0 + (en - st0 + 1) / 4 * 4, i, mm[4], p; \
	if (t < end) { \
		__m128i o_ = _mm_setr_epi32(0, m, 2 * m, 3 * m), m_ = _mm_set1_epi32(m), m4_ = _mm_set1_epi32(4 * m), mx_ = _mm_set1_epi32(KSW_NEG_INF), c_; \
		__m128i A_ = _mm_add_epi32(_mm_set1_epi32(m * (a + t)), o_), B_ = _mm_sub_epi32(_mm_set1_epi32(m * (b - t)), o_); \
		for (; t < end; t += 4) { /* the cells of r and r-1 at t..t+3 */ \
			__m128i h_ = KSW_EXTD2_BOUND_LOAD(H, t), p_; \
			p_ = _mm_add_epi32(h_, KSW_EXTD2_MIN32(A_, B_, c_)); \
			mx_ = KSW_EXTD2_MAX32(mx_, p_, c_); \
			p_ = _mm_add_epi32(_mm_sub_epi32(h_, ksw_i8x4_to_i32x4(&v8[t])), KSW_EXTD2_MIN32(_mm_add_epi32(A_, m_), B_, c_)); \
			mx_ = KSW_EXTD2_MAX32(mx_, p_, c_); \
			A_ = _mm_add_epi32(A_, m4_), B_ = _mm_sub_epi32(B_, m4_); \
		} \
		_mm_storeu_si128((__m128i*)mm, mx_); \
		for (i = 0; i < 4; ++i) M = M > mm[i]? M : mm[i]; \
	} \
	for (; t <= en; ++t) { \
		p = (H)[t] + m * (a + t < b - t? a + t : b - t), M = M > p? M : p; \
		p = (H)[t] - v8[t] + m * (a + 1 + t < b - t? a + 1 + t : b - t), M = M > p? M : p; \
	} \
	if (en < en0) p = (H)[en0] + m * (a + en0 < b - en0? a + en0 : b - en0), M = M > p? M : p; /* en0 of r only */ \
	if (lo1 < st0) p = (H)[lo1] + m * (a + 1 + lo1 < b - lo1? a + 1 + lo1 : b - lo1), M = M > p? M : p; /* st0-1 of r-1 only */ \
	return M;

static inline int32_t ksw_extd2_bound32(const int32_t *H, const int8_t *v8, int32_t st0, int32_t en0, int32_t lo1, int32_t hi1, int32_t a, int32_t b, int32_t m)
{
#define KSW_EXTD2_BOUND_LOAD(H, t) _mm_loadu_si128((const __m128i*)&(H)[t])
	KSW_EXTD2_BOUND_BODY(H)
#undef KSW_EXTD2_BOUND_LOAD
}

#ifdef KSW_H16
static inline int32_t ksw_extd2_bound16w(const int16_t *H, const int8_t *v8, int32_t st0, int32_t en0, int32_t lo1, int32_t hi1, int32_t a, int32_t b, int32_t m)
{
#define KSW_EXTD2_BOUND_LOAD(H, t) _mm_srai_epi32(_mm_unpacklo_epi16(_mm_loadl_epi64((const __m128i*)&(H)[t]), _mm_loadl_epi64((const __m128i*)&(H)[t])), 16)
	KSW_EXTD2_BOUND_BODY(H)
#undef KSW_EXTD2_BOUND_LOAD
}

/* The same over the int16 H in 8 int16 lanes, for m * (a + b + 2) < 32000 (else in 32-bit lanes): the matches-left terms then
 * fit, and H + that term saturates at 32767, which only raises the result; on the int16 path max <= 32767, so a saturated
 * result never stops the loop. */
static inline int32_t ksw_extd2_bound16(const int16_t *H, const int8_t *v8, int32_t st0, int32_t en0, int32_t lo1, int32_t hi1, int32_t a, int32_t b, int32_t m)
{
	int32_t t = st0, M = KSW_NEG_INF, en = hi1 < en0? hi1 : en0, end = st0 + (en - st0 + 1) / 8 * 8, p;
	if (m * (a + b + 2) >= 32000) return ksw_extd2_bound16w(H, v8, st0, en0, lo1, hi1, a, b, m);
	if (t < end) {
		__m128i o_ = _mm_setr_epi16(0, m, 2 * m, 3 * m, 4 * m, 5 * m, 6 * m, 7 * m), m_ = _mm_set1_epi16(m), m8_ = _mm_set1_epi16(8 * m);
		__m128i A_ = _mm_add_epi16(_mm_set1_epi16(m * (a + t)), o_), B_ = _mm_sub_epi16(_mm_set1_epi16(m * (b - t)), o_), mx_ = _mm_set1_epi16(-0x8000);
		for (; t < end; t += 8) { // the cells of r and r-1 at t..t+7
			__m128i h_ = _mm_loadu_si128((const __m128i*)&H[t]);
			mx_ = _mm_max_epi16(mx_, _mm_adds_epi16(h_, _mm_min_epi16(A_, B_)));
			h_ = _mm_sub_epi16(h_, ksw_i8x8_to_i16x8(_mm_loadl_epi64((const __m128i*)&v8[t])));
			mx_ = _mm_max_epi16(mx_, _mm_adds_epi16(h_, _mm_min_epi16(_mm_add_epi16(A_, m_), B_)));
			A_ = _mm_add_epi16(A_, m8_), B_ = _mm_sub_epi16(B_, m8_);
		}
		mx_ = _mm_max_epi16(mx_, _mm_srli_si128(mx_, 8));
		mx_ = _mm_max_epi16(mx_, _mm_srli_si128(mx_, 4));
		mx_ = _mm_max_epi16(mx_, _mm_srli_si128(mx_, 2));
		M = (int16_t)_mm_cvtsi128_si32(mx_);
	}
	for (; t <= en; ++t) {
		p = H[t] + m * (a + t < b - t? a + t : b - t), M = M > p? M : p;
		p = H[t] - v8[t] + m * (a + 1 + t < b - t? a + 1 + t : b - t), M = M > p? M : p;
	}
	if (en < en0) p = H[en0] + m * (a + en0 < b - en0? a + en0 : b - en0), M = M > p? M : p; // en0 of r only
	if (lo1 < st0) p = H[lo1] + m * (a + 1 + lo1 < b - lo1? a + 1 + lo1 : b - lo1), M = M > p? M : p; // st0-1 of r-1 only
	return M;
}
#endif
#undef KSW_EXTD2_BOUND_BODY

/* Calls with km == NULL take their working buffers (mem, H16 or H, mem2, off) from a scratch kept per thread between calls
 * and grown to the largest size asked for, instead of allocating and freeing them on every call; a thread's scratch is
 * released when the thread exits, by the destructor of a pthread key. A buffer above KSW_SCR_MAX bytes, a call with a km,
 * and a thread whose scratch cannot be set up allocate on every call (KSW_SCR_MAX=0 forces that, for testing). A scratch buffer
 * holds what an earlier call left in it: every buffer but mem is written before it is read, as with kmalloc, and mem's
 * regions after y2 are cleared, as kcalloc cleared them. */
#include <pthread.h>
#include <stdlib.h>

#ifndef KSW_SCR_MAX
#define KSW_SCR_MAX (4<<20)
#endif

enum { KSW_SCR_MEM, KSW_SCR_H16, KSW_SCR_H, KSW_SCR_MEM2, KSW_SCR_OFF, KSW_SCR_N };

typedef struct {
	void *raw[KSW_SCR_N], *buf[KSW_SCR_N]; // buf[i]: raw[i] aligned to 64 bytes
	size_t cap[KSW_SCR_N];
} ksw_scr_t;

static pthread_key_t ksw_scr_key;
static pthread_once_t ksw_scr_once = PTHREAD_ONCE_INIT;
static int ksw_scr_ok;

static void ksw_scr_destroy(void *x)
{
	ksw_scr_t *sc = (ksw_scr_t*)x;
	int i;
	for (i = 0; i < KSW_SCR_N; ++i) free(sc->raw[i]);
	free(sc);
}

static void ksw_scr_init(void)
{
	__atomic_store_n(&ksw_scr_ok, pthread_key_create(&ksw_scr_key, ksw_scr_destroy) == 0, __ATOMIC_RELEASE);
}

static __attribute__((noinline)) ksw_scr_t *ksw_scr_new(void) // this thread's first call: set up its scratch, or NULL
{
	ksw_scr_t *sc;
	if (pthread_once(&ksw_scr_once, ksw_scr_init) != 0 || !ksw_scr_ok) return 0;
	if ((sc = (ksw_scr_t*)pthread_getspecific(ksw_scr_key)) != 0) return sc;
	sc = (ksw_scr_t*)calloc(1, sizeof(ksw_scr_t));
	if (sc && pthread_setspecific(ksw_scr_key, sc) != 0) free(sc), sc = 0;
	return sc;
}

static inline ksw_scr_t *ksw_scr_get(void) // this thread's scratch, or NULL
{
	ksw_scr_t *sc;
	if (KSW_SCR_MAX <= 0) return 0;
	if (__atomic_load_n(&ksw_scr_ok, __ATOMIC_ACQUIRE) && (sc = (ksw_scr_t*)pthread_getspecific(ksw_scr_key)) != 0) return sc;
	return ksw_scr_new();
}

static __attribute__((noinline)) void *ksw_scr_grow(ksw_scr_t *sc, int i, size_t n) // slot i to at least n bytes, or NULL
{
	size_t m = n + (n >> 2);
	if (m > (size_t)KSW_SCR_MAX) m = KSW_SCR_MAX;
	free(sc->raw[i]);
	sc->raw[i] = malloc(m + 63);
	sc->buf[i] = sc->raw[i]? (void*)(((size_t)sc->raw[i] + 63) >> 6 << 6) : 0;
	sc->cap[i] = sc->raw[i]? m : 0;
	return sc->buf[i];
}

static inline void *ksw_scr_alloc(ksw_scr_t *sc, int i, size_t n) // n bytes of slot i, 64-byte aligned, or NULL, in which case the caller uses kmalloc() or kcalloc()
{
	if (sc == 0 || n == 0 || n > (size_t)KSW_SCR_MAX) return 0;
	return n <= sc->cap[i]? sc->buf[i] : ksw_scr_grow(sc, i, n);
}

static inline int ksw_scr_has(const ksw_scr_t *sc, int i, const void *p) // whether p is slot i's buffer
{
	return sc && p && p == sc->buf[i];
}

/* On x86-64 with gcc or clang, the function body below is compiled twice: as ksw_extd2_sse, and as
 * ksw_extd2_avx2 for the AVX2 target alone, which ksw_extd2_sse calls when the CPU and OS support AVX2.
 * The build flags stay SSE, and the SSE code is the fallback. The AVX2 copy runs the score fill and the
 * DP core on two 16-cell blocks per 256-bit operation, pairing the blocks from st_ and leaving an odd
 * last block to the 128-bit code; it computes the same 16-cell blocks and writes the same score bytes
 * (16-byte runs from st0) as the SSE code, so every cell, stored byte and result is the same.
 * KSW_NO_AVX2 leaves the AVX2 copy out, for testing. */
#if defined(__x86_64__) && defined(__GNUC__) && !defined(__ARM_NEON) && !defined(KSW_NO_AVX2)
#define KSW_EXTD2_AVX2 1
#include <immintrin.h>
#include <cpuid.h>
#define KSW_AVX2_FN __attribute__((target("avx2")))

static int ksw_extd2_has_avx2(void) // AVX2 in the CPU, and YMM state enabled by the OS; checked once
{
	static int avx2 = -1;
	int r = __atomic_load_n(&avx2, __ATOMIC_RELAXED);
	if (r < 0) {
		unsigned a, b, c, d, lo, hi;
		r = 0;
		if (__get_cpuid_max(0, 0) >= 7) {
			__cpuid(1, a, b, c, d);
			if ((c >> 27 & 1) && (c >> 28 & 1)) { // OSXSAVE and AVX
				__asm__ __volatile__ ("xgetbv" : "=a"(lo), "=d"(hi) : "c"(0));
				if ((lo & 6) == 6) { // XMM and YMM state
					__cpuid_count(7, 0, a, b, c, d);
					r = b >> 5 & 1; // AVX2
				}
			}
		}
		__atomic_store_n(&avx2, r, __ATOMIC_RELAXED);
	}
	return r;
}

#define KSW_EXTD2_BODY 256
#include "ksw2_extd2_sse.c" // static ksw_extd2_avx2()
#undef KSW_EXTD2_BODY
#endif

#define KSW_EXTD2_BODY 128
#include "ksw2_extd2_sse.c" // ksw_extd2_sse()
#undef KSW_EXTD2_BODY

#else // KSW_EXTD2_BODY: the function, included above once per vector width

#if KSW_EXTD2_BODY == 256
static KSW_AVX2_FN void ksw_extd2_avx2(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
									   int8_t q, int8_t e, int8_t q2, int8_t e2, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez)
#else
void ksw_extd2_sse(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
				   int8_t q, int8_t e, int8_t q2, int8_t e2, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez)
#endif
{
#if KSW_EXTD2_BODY == 128 && defined(KSW_EXTD2_AVX2)
	if (ksw_extd2_has_avx2()) {
		ksw_extd2_avx2(km, qlen, query, tlen, target, m, mat, q, e, q2, e2, w, zdrop, end_bonus, flag, ez);
		return;
	}
#endif
#define __dp_code_block1 \
	z = _mm_load_si128(&s[t]); \
	xt1 = _mm_load_si128(&x[t]);                     /* xt1 <- x[r-1][t..t+15] */ \
	tmp = _mm_srli_si128(xt1, 15);                   /* tmp <- x[r-1][t+15] */ \
	xt1 = _mm_or_si128(_mm_slli_si128(xt1, 1), x1_); /* xt1 <- x[r-1][t-1..t+14] */ \
	x1_ = tmp; \
	vt1 = _mm_load_si128(&v[t]);                     /* vt1 <- v[r-1][t..t+15] */ \
	tmp = _mm_srli_si128(vt1, 15);                   /* tmp <- v[r-1][t+15] */ \
	vt1 = _mm_or_si128(_mm_slli_si128(vt1, 1), v1_); /* vt1 <- v[r-1][t-1..t+14] */ \
	v1_ = tmp; \
	a = _mm_add_epi8(xt1, vt1);                      /* a <- x[r-1][t-1..t+14] + v[r-1][t-1..t+14] */ \
	ut = _mm_load_si128(&u[t]);                      /* ut <- u[t..t+15] */ \
	b = _mm_add_epi8(_mm_load_si128(&y[t]), ut);     /* b <- y[r-1][t..t+15] + u[r-1][t..t+15] */ \
	x2t1= _mm_load_si128(&x2[t]); \
	tmp = _mm_srli_si128(x2t1, 15); \
	x2t1= _mm_or_si128(_mm_slli_si128(x2t1, 1), x21_); \
	x21_= tmp; \
	a2= _mm_add_epi8(x2t1, vt1); \
	b2= _mm_add_epi8(_mm_load_si128(&y2[t]), ut);

#define __dp_code_block2 \
	_mm_store_si128(&u[t], _mm_sub_epi8(z, vt1));    /* u[r][t..t+15] <- z - v[r-1][t-1..t+14] */ \
	_mm_store_si128(&v[t], _mm_sub_epi8(z, ut));     /* v[r][t..t+15] <- z - u[r-1][t..t+15] */ \
	tmp = _mm_sub_epi8(z, q_); \
	a = _mm_sub_epi8(a, tmp); \
	b = _mm_sub_epi8(b, tmp); \
	tmp = _mm_sub_epi8(z, q2_); \
	a2= _mm_sub_epi8(a2, tmp); \
	b2= _mm_sub_epi8(b2, tmp);

#if KSW_EXTD2_BODY == 256
	// blocks t and t+1 at once; yx1_, yv1_ and yx21_ hold the last 32 bytes loaded, whose byte 31 is shifted in
#define __dp_code_block1_256 \
	z = _mm256_loadu_si256((__m256i*)&s[t]); \
	tmp = _mm256_loadu_si256((__m256i*)&x[t]);                                   /* tmp <- x[r-1][t..t+31] */ \
	xt1 = _mm256_alignr_epi8(tmp, _mm256_permute2x128_si256(tmp, yx1_, 0x03), 15); /* xt1 <- x[r-1][t-1..t+30] */ \
	yx1_ = tmp; \
	tmp = _mm256_loadu_si256((__m256i*)&v[t]); \
	vt1 = _mm256_alignr_epi8(tmp, _mm256_permute2x128_si256(tmp, yv1_, 0x03), 15); \
	yv1_ = tmp; \
	a = _mm256_add_epi8(xt1, vt1); \
	ut = _mm256_loadu_si256((__m256i*)&u[t]); \
	b = _mm256_add_epi8(_mm256_loadu_si256((__m256i*)&y[t]), ut); \
	tmp = _mm256_loadu_si256((__m256i*)&x2[t]); \
	x2t1 = _mm256_alignr_epi8(tmp, _mm256_permute2x128_si256(tmp, yx21_, 0x03), 15); \
	yx21_ = tmp; \
	a2 = _mm256_add_epi8(x2t1, vt1); \
	b2 = _mm256_add_epi8(_mm256_loadu_si256((__m256i*)&y2[t]), ut);

#define __dp_code_block2_256 \
	_mm256_storeu_si256((__m256i*)&u[t], _mm256_sub_epi8(z, vt1)); \
	_mm256_storeu_si256((__m256i*)&v[t], _mm256_sub_epi8(z, ut)); \
	tmp = _mm256_sub_epi8(z, yq_); \
	a = _mm256_sub_epi8(a, tmp); \
	b = _mm256_sub_epi8(b, tmp); \
	tmp = _mm256_sub_epi8(z, yq2_); \
	a2 = _mm256_sub_epi8(a2, tmp); \
	b2 = _mm256_sub_epi8(b2, tmp);
#endif

	int r, t, qe = q + e, n_col_, *off = 0, *off_end = 0, tlen_, qlen_, last_st, last_en, last_max_H = 0, wl, wr, max_sc, min_sc, long_thres, long_diff;
	int with_cigar = !(flag&KSW_EZ_SCORE_ONLY), approx_max = !!(flag&KSW_EZ_APPROX_MAX), zd_all, bound, bnext = 0;
	int32_t eb0 = end_bonus > 0? end_bonus : 0, bmax = -1;
	int32_t *H = 0, H0 = 0, last_H0_t = 0;
#ifdef KSW_H16
	int16_t *H16 = 0;
#endif
	uint8_t *qr, *sf, *mem, *mem2 = 0;
	ksw_scr_t *sc;
	size_t mem_n;
	__m128i q_, q2_, qe_, qe2_, zero_, sc_mch_, sc_mis_, m1_, sc_N_;
	__m128i *u, *v, *x, *y, *x2, *y2, *s, *p = 0;
#if KSW_EXTD2_BODY == 256
	__m256i yq_, yq2_, yqe_, yqe2_, yzero_, ysc_mch_, ysc_mis_, ym1_, ysc_N_;
#endif

	ksw_reset_extz(ez);
	if (m <= 1 || qlen <= 0 || tlen <= 0) return;

	if (q2 + e2 < q + e) t = q, q = q2, q2 = t, t = e, e = e2, e2 = t; // make sure q+e no larger than q2+e2

	zero_   = _mm_set1_epi8(0);
	q_      = _mm_set1_epi8(q);
	q2_     = _mm_set1_epi8(q2);
	qe_     = _mm_set1_epi8(q + e);
	qe2_    = _mm_set1_epi8(q2 + e2);
	sc_mch_ = _mm_set1_epi8(mat[0]);
	sc_mis_ = _mm_set1_epi8(mat[1]);
	sc_N_   = mat[m*m-1] == 0? _mm_set1_epi8(-e2) : _mm_set1_epi8(mat[m*m-1]);
	m1_     = _mm_set1_epi8(m - 1); // wildcard
#if KSW_EXTD2_BODY == 256
	yzero_   = _mm256_set1_epi8(0);
	yq_      = _mm256_set1_epi8(q);
	yq2_     = _mm256_set1_epi8(q2);
	yqe_     = _mm256_set1_epi8(q + e);
	yqe2_    = _mm256_set1_epi8(q2 + e2);
	ysc_mch_ = _mm256_set1_epi8(mat[0]);
	ysc_mis_ = _mm256_set1_epi8(mat[1]);
	ysc_N_   = mat[m*m-1] == 0? _mm256_set1_epi8(-e2) : _mm256_set1_epi8(mat[m*m-1]);
	ym1_     = _mm256_set1_epi8(m - 1);
#endif

	if (w < 0) w = tlen > qlen? tlen : qlen;
	wl = wr = w;
	tlen_ = (tlen + 15) / 16;
	n_col_ = qlen < tlen? qlen : tlen;
	n_col_ = ((n_col_ < w + 1? n_col_ : w + 1) + 15) / 16 + 1;
	qlen_ = (qlen + 15) / 16;
	for (t = 1, max_sc = mat[0], min_sc = mat[1]; t < m * m; ++t) {
		max_sc = max_sc > mat[t]? max_sc : mat[t];
		min_sc = min_sc < mat[t]? min_sc : mat[t];
	}
	if (-min_sc > 2 * (q + e)) return; // otherwise, we won't see any mismatches

	long_thres = e != e2? (q2 - q) / (e - e2) - 1 : 0;
	if (q2 + e2 + long_thres * e2 > q + e + long_thres * e)
		++long_thres;
	long_diff = long_thres * (e - e2) - (q2 - q) - e2;
	// ksw_apply_zdrop reads max_t when max_H > ez->max and otherwise only in a Z-drop test, ez->max - max_H > zdrop + l * e2
	// with l <= max(qlen, tlen), that cannot pass unless ez->max - max_H > zdrop, given e2 >= 0 and no overflow; else max_t always
	zd_all = zdrop >= 0 && (e2 < 0 || zdrop > INT32_MAX - (int64_t)(qlen > tlen? qlen : tlen) * e2);
#ifdef KSW_ALL_MAX_T
	zd_all = 1; // for testing: max_t on every anti-diagonal
#endif
	// KSW_EZ_EXTZ_BOUND, only where the argument at its stop below holds: gap costs >= 0 and a positive largest score, so that
	// the boundary row and column fall and the differences stay well inside int8; w >= 2 and qlen < tlen + w - 1, so that no
	// anti-diagonal's band is a single cell at offset t-i = 1-w, the one kind of cell whose H the exact-max block takes from
	// a stale H[en0-1]; and lengths for which the bound's int32 sums cannot overflow
	bound = (flag & KSW_EZ_EXTZ_BOUND) && (flag & KSW_EZ_EXTZ_ONLY) && !approx_max && max_sc > 0 && q >= 0 && e >= 0 && q2 >= 0 && e2 >= 0
		&& max_sc + q2 + e2 <= 63 && w >= 2 && (int64_t)tlen + w - 1 > qlen && qlen + tlen < 1<<22;

	sc = km == 0? ksw_scr_get() : 0;
	mem_n = (size_t)(tlen_ * 8 + qlen_ + 1) * 16;
	if ((mem = (uint8_t*)ksw_scr_alloc(sc, KSW_SCR_MEM, mem_n)) == 0)
		mem = (uint8_t*)kcalloc(km, tlen_ * 8 + qlen_ + 1, 16);
	u = (__m128i*)(((size_t)mem + 15) >> 4 << 4); // 16-byte aligned
	v = u + tlen_, x = v + tlen_, y = x + tlen_, x2 = y + tlen_, y2 = x2 + tlen_;
	s = y2 + tlen_, sf = (uint8_t*)(s + tlen_), qr = sf + tlen_ * 16;
	if (ksw_scr_has(sc, KSW_SCR_MEM, mem)) memset(s, 0, mem + mem_n - (uint8_t*)s); // as kcalloc left them; u..y2 are set below
	memset(u,  -q  - e,  tlen_ * 16);
	memset(v,  -q  - e,  tlen_ * 16);
	memset(x,  -q  - e,  tlen_ * 16);
	memset(y,  -q  - e,  tlen_ * 16);
	memset(x2, -q2 - e2, tlen_ * 16);
	memset(y2, -q2 - e2, tlen_ * 16);
	if (!approx_max) {
#ifdef KSW_H16
		if (tlen < KSW_H16_MAX_T) { // no initialization: only entries 0..en0 are read, and each was written when it was en0
			if ((H16 = (int16_t*)ksw_scr_alloc(sc, KSW_SCR_H16, tlen_ * 16 * 2)) == 0)
				H16 = (int16_t*)kmalloc(km, tlen_ * 16 * 2);
		} else
#endif
		{
			if ((H = (int32_t*)ksw_scr_alloc(sc, KSW_SCR_H, tlen_ * 16 * 4)) == 0)
				H = (int32_t*)kmalloc(km, tlen_ * 16 * 4);
			for (t = 0; t < tlen_ * 16; ++t) H[t] = KSW_NEG_INF;
		}
	}
	if (with_cigar) { // mem2 and off: no initialization, the backtrack reads only cells and rows the DP wrote
		if ((mem2 = (uint8_t*)ksw_scr_alloc(sc, KSW_SCR_MEM2, ((size_t)(qlen + tlen - 1) * n_col_ + 1) * 16)) == 0)
			mem2 = (uint8_t*)kmalloc(km, ((size_t)(qlen + tlen - 1) * n_col_ + 1) * 16);
		p = (__m128i*)(((size_t)mem2 + 15) >> 4 << 4);
		if ((off = (int*)ksw_scr_alloc(sc, KSW_SCR_OFF, (qlen + tlen - 1) * sizeof(int) * 2)) == 0)
			off = (int*)kmalloc(km, (qlen + tlen - 1) * sizeof(int) * 2);
		off_end = off + qlen + tlen - 1;
	}

	for (t = 0; t < qlen; ++t) qr[t] = query[qlen - 1 - t];
	memcpy(sf, target, tlen);

	for (r = 0, last_st = last_en = -1; r < qlen + tlen - 1; ++r) {
		int st = 0, en = tlen - 1, st0, en0, st_, en_;
		int8_t x1, x21, v1;
		uint8_t *qrr = qr + (qlen - 1 - r);
		int8_t *u8 = (int8_t*)u, *v8 = (int8_t*)v, *x8 = (int8_t*)x, *x28 = (int8_t*)x2;
		__m128i x1_, x21_, v1_;
		// find the boundaries
		if (st < r - qlen + 1) st = r - qlen + 1;
		if (en > r) en = r;
		if (st < (r-wr+1)>>1) st = (r-wr+1)>>1; // take the ceil
		if (en > (r+wl)>>1) en = (r+wl)>>1; // take the floor
		if (st > en) {
			ez->zdropped = 1;
			break;
		}
		st0 = st, en0 = en;
		st = st / 16 * 16, en = (en + 16) / 16 * 16 - 1;
		// set boundary conditions
		if (st > 0) {
			if (st - 1 >= last_st && st - 1 <= last_en) {
				x1 = x8[st - 1], x21 = x28[st - 1], v1 = v8[st - 1]; // (r-1,s-1) calculated in the last round
			} else {
				x1 = -q - e, x21 = -q2 - e2;
				v1 = -q - e;
			}
		} else {
			x1 = -q - e, x21 = -q2 - e2;
			v1 = r == 0? -q - e : r < long_thres? -e : r == long_thres? long_diff : -e2;
		}
		if (en >= r) {
			((int8_t*)y)[r] = -q - e, ((int8_t*)y2)[r] = -q2 - e2;
			u8[r] = r == 0? -q - e : r < long_thres? -e : r == long_thres? long_diff : -e2;
		}
		// loop fission: set scores first
		if (!(flag & KSW_EZ_GENERIC_SC)) {
			t = st0;
#if KSW_EXTD2_BODY == 256
			for (; t + 16 <= en0; t += 32) { // two of the 16-byte steps below at once
				__m256i sq, st, tmp, mask;
				sq = _mm256_loadu_si256((__m256i*)&sf[t]);
				st = _mm256_loadu_si256((__m256i*)&qrr[t]);
				mask = _mm256_or_si256(_mm256_cmpeq_epi8(sq, ym1_), _mm256_cmpeq_epi8(st, ym1_));
				tmp = _mm256_cmpeq_epi8(sq, st);
				tmp = _mm256_blendv_epi8(ysc_mis_, ysc_mch_, tmp);
				tmp = _mm256_blendv_epi8(tmp,      ysc_N_,   mask);
				_mm256_storeu_si256((__m256i*)((int8_t*)s + t), tmp);
			}
#endif
			for (; t <= en0; t += 16) {
				__m128i sq, st, tmp, mask;
				sq = _mm_loadu_si128((__m128i*)&sf[t]);
				st = _mm_loadu_si128((__m128i*)&qrr[t]);
				mask = _mm_or_si128(_mm_cmpeq_epi8(sq, m1_), _mm_cmpeq_epi8(st, m1_));
				tmp = _mm_cmpeq_epi8(sq, st);
#ifdef __SSE4_1__
				tmp = _mm_blendv_epi8(sc_mis_, sc_mch_, tmp);
				tmp = _mm_blendv_epi8(tmp,     sc_N_,   mask);
#else
				tmp = _mm_or_si128(_mm_andnot_si128(tmp,  sc_mis_), _mm_and_si128(tmp,  sc_mch_));
				tmp = _mm_or_si128(_mm_andnot_si128(mask, tmp),     _mm_and_si128(mask, sc_N_));
#endif
				_mm_storeu_si128((__m128i*)((int8_t*)s + t), tmp);
			}
		} else {
			for (t = st0; t <= en0; ++t)
				((uint8_t*)s)[t] = mat[sf[t] * m + qrr[t]];
		}
		// core loop
		x1_  = _mm_cvtsi32_si128((uint8_t)x1);
		x21_ = _mm_cvtsi32_si128((uint8_t)x21);
		v1_  = _mm_cvtsi32_si128((uint8_t)v1);
		st_ = st / 16, en_ = en / 16;
		assert(en_ - st_ + 1 <= n_col_);
		t = st_;
#if KSW_EXTD2_BODY == 256
		if (t < en_) { // blocks in pairs from st_; an odd last block goes to the 128-bit loops below
			__m256i yx1_ = _mm256_set1_epi8(x1), yx21_ = _mm256_set1_epi8(x21), yv1_ = _mm256_set1_epi8(v1); // byte 31 is shifted in
			if (!with_cigar) { // score only
				for (; t < en_; t += 2) {
					__m256i z, a, b, a2, b2, xt1, x2t1, vt1, ut, tmp;
					__dp_code_block1_256;
					z = _mm256_max_epi8(z, a);
					z = _mm256_max_epi8(z, b);
					z = _mm256_max_epi8(z, a2);
					z = _mm256_max_epi8(z, b2);
					z = _mm256_min_epi8(z, ysc_mch_);
					__dp_code_block2_256;
					_mm256_storeu_si256((__m256i*)&x[t],  _mm256_sub_epi8(_mm256_max_epi8(a,  yzero_), yqe_));
					_mm256_storeu_si256((__m256i*)&y[t],  _mm256_sub_epi8(_mm256_max_epi8(b,  yzero_), yqe_));
					_mm256_storeu_si256((__m256i*)&x2[t], _mm256_sub_epi8(_mm256_max_epi8(a2, yzero_), yqe2_));
					_mm256_storeu_si256((__m256i*)&y2[t], _mm256_sub_epi8(_mm256_max_epi8(b2, yzero_), yqe2_));
				}
			} else if (!(flag&KSW_EZ_RIGHT)) { // gap left-alignment
				__m128i *pr = p + (size_t)r * n_col_ - st_;
				for (; t < en_; t += 2) {
					__m256i d, z, a, b, a2, b2, xt1, x2t1, vt1, ut, tmp;
					__dp_code_block1_256;
					d = _mm256_and_si256(_mm256_cmpgt_epi8(a, z), _mm256_set1_epi8(1));       // d = a  > z? 1 : 0
					z = _mm256_max_epi8(z, a);
					d = _mm256_blendv_epi8(d, _mm256_set1_epi8(2), _mm256_cmpgt_epi8(b,  z)); // d = b  > z? 2 : d
					z = _mm256_max_epi8(z, b);
					d = _mm256_blendv_epi8(d, _mm256_set1_epi8(3), _mm256_cmpgt_epi8(a2, z)); // d = a2 > z? 3 : d
					z = _mm256_max_epi8(z, a2);
					d = _mm256_blendv_epi8(d, _mm256_set1_epi8(4), _mm256_cmpgt_epi8(b2, z)); // d = b2 > z? 4 : d
					z = _mm256_max_epi8(z, b2);
					z = _mm256_min_epi8(z, ysc_mch_);
					__dp_code_block2_256;
					tmp = _mm256_cmpgt_epi8(a, yzero_);
					_mm256_storeu_si256((__m256i*)&x[t],  _mm256_sub_epi8(_mm256_and_si256(tmp, a),  yqe_));
					d = _mm256_or_si256(d, _mm256_and_si256(tmp, _mm256_set1_epi8(0x08)));
					tmp = _mm256_cmpgt_epi8(b, yzero_);
					_mm256_storeu_si256((__m256i*)&y[t],  _mm256_sub_epi8(_mm256_and_si256(tmp, b),  yqe_));
					d = _mm256_or_si256(d, _mm256_and_si256(tmp, _mm256_set1_epi8(0x10)));
					tmp = _mm256_cmpgt_epi8(a2, yzero_);
					_mm256_storeu_si256((__m256i*)&x2[t], _mm256_sub_epi8(_mm256_and_si256(tmp, a2), yqe2_));
					d = _mm256_or_si256(d, _mm256_and_si256(tmp, _mm256_set1_epi8(0x20)));
					tmp = _mm256_cmpgt_epi8(b2, yzero_);
					_mm256_storeu_si256((__m256i*)&y2[t], _mm256_sub_epi8(_mm256_and_si256(tmp, b2), yqe2_));
					d = _mm256_or_si256(d, _mm256_and_si256(tmp, _mm256_set1_epi8(0x40)));
					_mm256_storeu_si256((__m256i*)&pr[t], d);
				}
			} else { // gap right-alignment
				__m128i *pr = p + (size_t)r * n_col_ - st_;
				for (; t < en_; t += 2) {
					__m256i d, z, a, b, a2, b2, xt1, x2t1, vt1, ut, tmp;
					__dp_code_block1_256;
					d = _mm256_andnot_si256(_mm256_cmpgt_epi8(z, a), _mm256_set1_epi8(1));    // d = z > a?  0 : 1
					z = _mm256_max_epi8(z, a);
					d = _mm256_blendv_epi8(_mm256_set1_epi8(2), d, _mm256_cmpgt_epi8(z, b));  // d = z > b?  d : 2
					z = _mm256_max_epi8(z, b);
					d = _mm256_blendv_epi8(_mm256_set1_epi8(3), d, _mm256_cmpgt_epi8(z, a2)); // d = z > a2? d : 3
					z = _mm256_max_epi8(z, a2);
					d = _mm256_blendv_epi8(_mm256_set1_epi8(4), d, _mm256_cmpgt_epi8(z, b2)); // d = z > b2? d : 4
					z = _mm256_max_epi8(z, b2);
					z = _mm256_min_epi8(z, ysc_mch_);
					__dp_code_block2_256;
					tmp = _mm256_cmpgt_epi8(yzero_, a);
					_mm256_storeu_si256((__m256i*)&x[t],  _mm256_sub_epi8(_mm256_andnot_si256(tmp, a),  yqe_));
					d = _mm256_or_si256(d, _mm256_andnot_si256(tmp, _mm256_set1_epi8(0x08)));
					tmp = _mm256_cmpgt_epi8(yzero_, b);
					_mm256_storeu_si256((__m256i*)&y[t],  _mm256_sub_epi8(_mm256_andnot_si256(tmp, b),  yqe_));
					d = _mm256_or_si256(d, _mm256_andnot_si256(tmp, _mm256_set1_epi8(0x10)));
					tmp = _mm256_cmpgt_epi8(yzero_, a2);
					_mm256_storeu_si256((__m256i*)&x2[t], _mm256_sub_epi8(_mm256_andnot_si256(tmp, a2), yqe2_));
					d = _mm256_or_si256(d, _mm256_andnot_si256(tmp, _mm256_set1_epi8(0x20)));
					tmp = _mm256_cmpgt_epi8(yzero_, b2);
					_mm256_storeu_si256((__m256i*)&y2[t], _mm256_sub_epi8(_mm256_andnot_si256(tmp, b2), yqe2_));
					d = _mm256_or_si256(d, _mm256_andnot_si256(tmp, _mm256_set1_epi8(0x40)));
					_mm256_storeu_si256((__m256i*)&pr[t], d);
				}
			}
			x1_  = _mm_srli_si128(_mm256_extracti128_si256(yx1_,  1), 15); // the 128-bit carries: byte 31 of the last pair loaded
			x21_ = _mm_srli_si128(_mm256_extracti128_si256(yx21_, 1), 15);
			v1_  = _mm_srli_si128(_mm256_extracti128_si256(yv1_,  1), 15);
		}
#endif
		if (!with_cigar) { // score only
			for (; t <= en_; ++t) {
				__m128i z, a, b, a2, b2, xt1, x2t1, vt1, ut, tmp;
				__dp_code_block1;
#ifdef __SSE4_1__
				z = _mm_max_epi8(z, a);
				z = _mm_max_epi8(z, b);
				z = _mm_max_epi8(z, a2);
				z = _mm_max_epi8(z, b2);
				z = _mm_min_epi8(z, sc_mch_);
				__dp_code_block2; // save u[] and v[]; update a, b, a2 and b2
				_mm_store_si128(&x[t],  _mm_sub_epi8(_mm_max_epi8(a,  zero_), qe_));
				_mm_store_si128(&y[t],  _mm_sub_epi8(_mm_max_epi8(b,  zero_), qe_));
				_mm_store_si128(&x2[t], _mm_sub_epi8(_mm_max_epi8(a2, zero_), qe2_));
				_mm_store_si128(&y2[t], _mm_sub_epi8(_mm_max_epi8(b2, zero_), qe2_));
#else
				tmp = _mm_cmpgt_epi8(a,  z);
				z = _mm_or_si128(_mm_andnot_si128(tmp, z), _mm_and_si128(tmp, a));
				tmp = _mm_cmpgt_epi8(b,  z);
				z = _mm_or_si128(_mm_andnot_si128(tmp, z), _mm_and_si128(tmp, b));
				tmp = _mm_cmpgt_epi8(a2, z);
				z = _mm_or_si128(_mm_andnot_si128(tmp, z), _mm_and_si128(tmp, a2));
				tmp = _mm_cmpgt_epi8(b2, z);
				z = _mm_or_si128(_mm_andnot_si128(tmp, z), _mm_and_si128(tmp, b2));
				tmp = _mm_cmplt_epi8(sc_mch_, z);
				z = _mm_or_si128(_mm_and_si128(tmp, sc_mch_), _mm_andnot_si128(tmp, z));
				__dp_code_block2;
				tmp = _mm_cmpgt_epi8(a, zero_);
				_mm_store_si128(&x[t],  _mm_sub_epi8(_mm_and_si128(tmp, a),  qe_));
				tmp = _mm_cmpgt_epi8(b, zero_);
				_mm_store_si128(&y[t],  _mm_sub_epi8(_mm_and_si128(tmp, b),  qe_));
				tmp = _mm_cmpgt_epi8(a2, zero_);
				_mm_store_si128(&x2[t], _mm_sub_epi8(_mm_and_si128(tmp, a2), qe2_));
				tmp = _mm_cmpgt_epi8(b2, zero_);
				_mm_store_si128(&y2[t], _mm_sub_epi8(_mm_and_si128(tmp, b2), qe2_));
#endif
			}
		} else if (!(flag&KSW_EZ_RIGHT)) { // gap left-alignment
			__m128i *pr = p + (size_t)r * n_col_ - st_;
			off[r] = st, off_end[r] = en;
			for (; t <= en_; ++t) {
				__m128i d, z, a, b, a2, b2, xt1, x2t1, vt1, ut, tmp;
				__dp_code_block1;
#ifdef __SSE4_1__
				d = _mm_and_si128(_mm_cmpgt_epi8(a, z), _mm_set1_epi8(1));       // d = a  > z? 1 : 0
				z = _mm_max_epi8(z, a);
				d = _mm_blendv_epi8(d, _mm_set1_epi8(2), _mm_cmpgt_epi8(b,  z)); // d = b  > z? 2 : d
				z = _mm_max_epi8(z, b);
				d = _mm_blendv_epi8(d, _mm_set1_epi8(3), _mm_cmpgt_epi8(a2, z)); // d = a2 > z? 3 : d
				z = _mm_max_epi8(z, a2);
				d = _mm_blendv_epi8(d, _mm_set1_epi8(4), _mm_cmpgt_epi8(b2, z)); // d = a2 > z? 3 : d
				z = _mm_max_epi8(z, b2);
				z = _mm_min_epi8(z, sc_mch_);
#else // we need to emulate SSE4.1 intrinsics _mm_max_epi8() and _mm_blendv_epi8()
				tmp = _mm_cmpgt_epi8(a,  z);
				d = _mm_and_si128(tmp, _mm_set1_epi8(1));
				z = _mm_or_si128(_mm_andnot_si128(tmp, z), _mm_and_si128(tmp, a));
				tmp = _mm_cmpgt_epi8(b,  z);
				d = _mm_or_si128(_mm_andnot_si128(tmp, d), _mm_and_si128(tmp, _mm_set1_epi8(2)));
				z = _mm_or_si128(_mm_andnot_si128(tmp, z), _mm_and_si128(tmp, b));
				tmp = _mm_cmpgt_epi8(a2, z);
				d = _mm_or_si128(_mm_andnot_si128(tmp, d), _mm_and_si128(tmp, _mm_set1_epi8(3)));
				z = _mm_or_si128(_mm_andnot_si128(tmp, z), _mm_and_si128(tmp, a2));
				tmp = _mm_cmpgt_epi8(b2, z);
				d = _mm_or_si128(_mm_andnot_si128(tmp, d), _mm_and_si128(tmp, _mm_set1_epi8(4)));
				z = _mm_or_si128(_mm_andnot_si128(tmp, z), _mm_and_si128(tmp, b2));
				tmp = _mm_cmplt_epi8(sc_mch_, z);
				z = _mm_or_si128(_mm_and_si128(tmp, sc_mch_), _mm_andnot_si128(tmp, z));
#endif
				__dp_code_block2;
				tmp = _mm_cmpgt_epi8(a, zero_);
				_mm_store_si128(&x[t],  _mm_sub_epi8(_mm_and_si128(tmp, a),  qe_));
				d = _mm_or_si128(d, _mm_and_si128(tmp, _mm_set1_epi8(0x08))); // d = a > 0? 1<<3 : 0
				tmp = _mm_cmpgt_epi8(b, zero_);
				_mm_store_si128(&y[t],  _mm_sub_epi8(_mm_and_si128(tmp, b),  qe_));
				d = _mm_or_si128(d, _mm_and_si128(tmp, _mm_set1_epi8(0x10))); // d = b > 0? 1<<4 : 0
				tmp = _mm_cmpgt_epi8(a2, zero_);
				_mm_store_si128(&x2[t], _mm_sub_epi8(_mm_and_si128(tmp, a2), qe2_));
				d = _mm_or_si128(d, _mm_and_si128(tmp, _mm_set1_epi8(0x20))); // d = a > 0? 1<<5 : 0
				tmp = _mm_cmpgt_epi8(b2, zero_);
				_mm_store_si128(&y2[t], _mm_sub_epi8(_mm_and_si128(tmp, b2), qe2_));
				d = _mm_or_si128(d, _mm_and_si128(tmp, _mm_set1_epi8(0x40))); // d = b > 0? 1<<6 : 0
				_mm_store_si128(&pr[t], d);
			}
		} else { // gap right-alignment
			__m128i *pr = p + (size_t)r * n_col_ - st_;
			off[r] = st, off_end[r] = en;
			for (; t <= en_; ++t) {
				__m128i d, z, a, b, a2, b2, xt1, x2t1, vt1, ut, tmp;
				__dp_code_block1;
#ifdef __SSE4_1__
				d = _mm_andnot_si128(_mm_cmpgt_epi8(z, a), _mm_set1_epi8(1));    // d = z > a?  0 : 1
				z = _mm_max_epi8(z, a);
				d = _mm_blendv_epi8(_mm_set1_epi8(2), d, _mm_cmpgt_epi8(z, b));  // d = z > b?  d : 2
				z = _mm_max_epi8(z, b);
				d = _mm_blendv_epi8(_mm_set1_epi8(3), d, _mm_cmpgt_epi8(z, a2)); // d = z > a2? d : 3
				z = _mm_max_epi8(z, a2);
				d = _mm_blendv_epi8(_mm_set1_epi8(4), d, _mm_cmpgt_epi8(z, b2)); // d = z > b2? d : 4
				z = _mm_max_epi8(z, b2);
				z = _mm_min_epi8(z, sc_mch_);
#else // we need to emulate SSE4.1 intrinsics _mm_max_epi8() and _mm_blendv_epi8()
				tmp = _mm_cmpgt_epi8(z, a);
				d = _mm_andnot_si128(tmp, _mm_set1_epi8(1));
				z = _mm_or_si128(_mm_and_si128(tmp, z), _mm_andnot_si128(tmp, a));
				tmp = _mm_cmpgt_epi8(z, b);
				d = _mm_or_si128(_mm_and_si128(tmp, d), _mm_andnot_si128(tmp, _mm_set1_epi8(2)));
				z = _mm_or_si128(_mm_and_si128(tmp, z), _mm_andnot_si128(tmp, b));
				tmp = _mm_cmpgt_epi8(z, a2);
				d = _mm_or_si128(_mm_and_si128(tmp, d), _mm_andnot_si128(tmp, _mm_set1_epi8(3)));
				z = _mm_or_si128(_mm_and_si128(tmp, z), _mm_andnot_si128(tmp, a2));
				tmp = _mm_cmpgt_epi8(z, b2);
				d = _mm_or_si128(_mm_and_si128(tmp, d), _mm_andnot_si128(tmp, _mm_set1_epi8(4)));
				z = _mm_or_si128(_mm_and_si128(tmp, z), _mm_andnot_si128(tmp, b2));
				tmp = _mm_cmplt_epi8(sc_mch_, z);
				z = _mm_or_si128(_mm_and_si128(tmp, sc_mch_), _mm_andnot_si128(tmp, z));
#endif
				__dp_code_block2;
				tmp = _mm_cmpgt_epi8(zero_, a);
				_mm_store_si128(&x[t],  _mm_sub_epi8(_mm_andnot_si128(tmp, a),  qe_));
				d = _mm_or_si128(d, _mm_andnot_si128(tmp, _mm_set1_epi8(0x08))); // d = a > 0? 1<<3 : 0
				tmp = _mm_cmpgt_epi8(zero_, b);
				_mm_store_si128(&y[t],  _mm_sub_epi8(_mm_andnot_si128(tmp, b),  qe_));
				d = _mm_or_si128(d, _mm_andnot_si128(tmp, _mm_set1_epi8(0x10))); // d = b > 0? 1<<4 : 0
				tmp = _mm_cmpgt_epi8(zero_, a2);
				_mm_store_si128(&x2[t], _mm_sub_epi8(_mm_andnot_si128(tmp, a2), qe2_));
				d = _mm_or_si128(d, _mm_andnot_si128(tmp, _mm_set1_epi8(0x20))); // d = a > 0? 1<<5 : 0
				tmp = _mm_cmpgt_epi8(zero_, b2);
				_mm_store_si128(&y2[t], _mm_sub_epi8(_mm_andnot_si128(tmp, b2), qe2_));
				d = _mm_or_si128(d, _mm_andnot_si128(tmp, _mm_set1_epi8(0x40))); // d = b > 0? 1<<6 : 0
				_mm_store_si128(&pr[t], d);
			}
		}
		if (!approx_max) { // find the exact max with a 32-bit score array (16-bit while the values fit; KSW_H16)
			// max_t only feeds ksw_apply_zdrop, so it is found only on the anti-diagonals where that reads it (need_t): the
			// first pass gives max_H, max_t is -1 unless a scalar tail cell raises the maximum, and the lanes' first maxima
			// are found by a second pass over the stored H, the first pass's tie-break code unchanged
			int32_t max_H, max_t, h_en0, h_st0, need_t;
#ifdef KSW_H16
			if (H16) {
				int32_t lo;
				__m128i min_H_ = _mm_setzero_si128();
				if (r > 0) {
					// the 32-bit block's result exactly: its cells st0..en1-1 in 4 lanes (4-lane j = int16 lanes j and j+4, as
					// (t-st0)%4 == (t-st0)%8 % 4), H[en0] seeding the max, each lane's first maximum, lanes 0-3 reduced with
					// strict >, then the same scalar tail
					int32_t en1 = st0 + (en0 - st0) / 4 * 4, en2 = st0 + (en0 - st0) / 8 * 8, M, j;
					__m128i max_H_, m_;
					max_H = H16[en0] = en0 > 0? H16[en0-1] + u8[en0] : H16[en0] + v8[en0]; // special casing the last element
					max_H_ = min_H_ = _mm_set1_epi16(max_H);
					for (t = st0; t < en2; t += 8) { // H[t]+=v8[t]; the lanes' maximum and minimum
						__m128i H1;
						H1 = _mm_loadu_si128((__m128i*)&H16[t]);
						H1 = _mm_add_epi16(H1, ksw_i8x8_to_i16x8(_mm_loadl_epi64((__m128i*)&v8[t])));
						_mm_storeu_si128((__m128i*)&H16[t], H1);
						max_H_ = _mm_max_epi16(max_H_, H1);
						min_H_ = _mm_min_epi16(min_H_, H1);
					}
					if (t < en1) { // four more cells, in lanes 0-3; lanes 4-7 hold 0 (harmless to the minimum) and then -32768 (never a maximum)
						__m128i H1;
						H1 = _mm_loadl_epi64((__m128i*)&H16[t]);
						H1 = _mm_add_epi16(H1, ksw_i8x8_to_i16x8(_mm_cvtsi32_si128(*(int32_t*)&v8[t])));
						_mm_storel_epi64((__m128i*)&H16[t], H1);
						min_H_ = _mm_min_epi16(min_H_, H1);
						H1 = _mm_or_si128(H1, _mm_setr_epi16(0, 0, 0, 0, -0x8000, -0x8000, -0x8000, -0x8000));
						max_H_ = _mm_max_epi16(max_H_, H1);
						t += 4;
					}
					// M: the lanes' maximum (at least H[en0])
					m_ = _mm_max_epi16(max_H_, _mm_srli_si128(max_H_, 8));
					m_ = _mm_max_epi16(m_, _mm_srli_si128(m_, 4));
					m_ = _mm_max_epi16(m_, _mm_srli_si128(m_, 2));
					M = (int16_t)_mm_cvtsi128_si32(m_);
					max_H = M, max_t = -1;
					lo = H16[en0];
					for (; t < en0; ++t) { // for the rest of values that haven't been computed with SSE
						int32_t h = H16[t] + v8[t];
						H16[t] = h;
						if (h > max_H)
							max_H = h, max_t = t;
						lo = lo < h? lo : h;
					}
					need_t = ksw_extd2_need_max_t(ez, max_H, zdrop, zd_all);
					if (need_t && max_t < 0) { // M is the maximum: j, the lowest 4-lane holding it, and tt[j], its first position in t
						int16_t tt[8];
						__m128i max_t_, t_, t8_ = _mm_set1_epi16(8), eq, tv;
						max_H_ = _mm_set1_epi16(H16[en0]);
						max_t_ = _mm_set1_epi16(en0);
						t_ = _mm_set1_epi16(st0);
						for (t = st0; t < en2; t += 8) { // per lane: if(H[t]>max_H) max_H=H[t],max_t=t
							__m128i H1, tmp;
							H1 = _mm_loadu_si128((__m128i*)&H16[t]);
							tmp = _mm_cmpgt_epi16(H1, max_H_);
							max_H_ = _mm_max_epi16(max_H_, H1);
							max_t_ = ksw_blend16(max_t_, t_, tmp);
							t_ = _mm_add_epi16(t_, t8_);
						}
						if (t < en1) {
							__m128i H1, tmp;
							H1 = _mm_or_si128(_mm_loadl_epi64((__m128i*)&H16[t]), _mm_setr_epi16(0, 0, 0, 0, -0x8000, -0x8000, -0x8000, -0x8000));
							tmp = _mm_cmpgt_epi16(H1, max_H_);
							max_H_ = _mm_max_epi16(max_H_, H1);
							max_t_ = ksw_blend16(max_t_, t_, tmp);
						}
						eq = _mm_cmpeq_epi16(max_H_, _mm_set1_epi16(M));
						tv = _mm_add_epi16(max_t_, _mm_setr_epi16(0, 1, 2, 3, 4, 5, 6, 7));
						tv = ksw_blend16(_mm_set1_epi16(0x7fff), tv, eq);
						tv = _mm_min_epi16(tv, _mm_srli_si128(tv, 8));
						_mm_storeu_si128((__m128i*)tt, tv);
						j = __builtin_ctz(_mm_movemask_epi8(_mm_or_si128(eq, _mm_srli_si128(eq, 8))) & 0xff) >> 1;
						max_t = tt[j]; // M == H[en0]: every lane still holds (H[en0], en0), so j == 0 and tt[0] == en0
					}
				} else {
					H16[0] = v8[0] - qe, max_H = H16[0], max_t = 0, lo = max_H; // special casing r==0
					need_t = ksw_extd2_need_max_t(ez, max_H, zdrop, zd_all);
				}
				if (max_H > KSW_H16_LIM || lo < -KSW_H16_LIM || _mm_movemask_epi8(_mm_cmpgt_epi16(_mm_set1_epi16(-KSW_H16_LIM), min_H_))) {
					// a value near the int16 limits: continue with the 32-bit array, holding what it would hold now
					if ((H = (int32_t*)ksw_scr_alloc(sc, KSW_SCR_H, tlen_ * 16 * 4)) == 0)
						H = (int32_t*)kmalloc(km, tlen_ * 16 * 4);
					for (t = 0; t <= en0; ++t) H[t] = H16[t];
					for (; t < tlen_ * 16; ++t) H[t] = KSW_NEG_INF;
					if (!ksw_scr_has(sc, KSW_SCR_H16, H16)) kfree(km, H16);
					H16 = 0;
				}
			} else
#endif
			if (r > 0) {
				int32_t HH[4], en1 = st0 + (en0 - st0) / 4 * 4, i;
				__m128i max_H_;
				max_H = H[en0] = en0 > 0? H[en0-1] + u8[en0] : H[en0] + v8[en0]; // special casing the last element
				max_H_ = _mm_set1_epi32(max_H);
				for (t = st0; t < en1; t += 4) { // this implements: H[t]+=v8[t]-qe; and the lanes' maximum
					__m128i H1, v_;
					H1 = _mm_loadu_si128((__m128i*)&H[t]);
					v_ = ksw_i8x4_to_i32x4(&v8[t]);
					H1 = _mm_add_epi32(H1, v_);
					_mm_storeu_si128((__m128i*)&H[t], H1);
#ifdef __SSE4_1__
					max_H_ = _mm_max_epi32(max_H_, H1);
#else
					v_ = _mm_cmpgt_epi32(H1, max_H_);
					max_H_ = _mm_or_si128(_mm_and_si128(v_, H1), _mm_andnot_si128(v_, max_H_));
#endif
				}
				_mm_storeu_si128((__m128i*)HH, max_H_);
				for (i = 0; i < 4; ++i)
					max_H = max_H > HH[i]? max_H : HH[i];
				max_t = -1;
				for (; t < en0; ++t) { // for the rest of values that haven't been computed with SSE
					H[t] += (int32_t)v8[t];
					if (H[t] > max_H)
						max_H = H[t], max_t = t;
				}
				need_t = ksw_extd2_need_max_t(ez, max_H, zdrop, zd_all);
				if (need_t && max_t < 0) { // the lanes' first maxima over the stored H, lanes 0-3 reduced with strict >
					int32_t tt[4], mh = H[en0];
					__m128i max_t_, t_, t4_ = _mm_set1_epi32(4);
					max_t = en0;
					max_H_ = _mm_set1_epi32(mh);
					max_t_ = _mm_set1_epi32(max_t);
					t_ = _mm_set1_epi32(st0);
					for (t = st0; t < en1; t += 4) { // per lane: if(H[t]>max_H) max_H=H[t],max_t=t
						__m128i H1, tmp;
						H1 = _mm_loadu_si128((__m128i*)&H[t]);
						tmp = _mm_cmpgt_epi32(H1, max_H_);
#ifdef __SSE4_1__
						max_H_ = _mm_max_epi32(max_H_, H1);
						max_t_ = _mm_blendv_epi8(max_t_, t_, tmp);
#else
						max_H_ = _mm_or_si128(_mm_and_si128(tmp, H1), _mm_andnot_si128(tmp, max_H_));
						max_t_ = _mm_or_si128(_mm_and_si128(tmp, t_), _mm_andnot_si128(tmp, max_t_));
#endif
						t_ = _mm_add_epi32(t_, t4_);
					}
					_mm_storeu_si128((__m128i*)HH, max_H_);
					_mm_storeu_si128((__m128i*)tt, max_t_);
					for (i = 0; i < 4; ++i)
						if (mh < HH[i]) mh = HH[i], max_t = tt[i] + i;
				}
			} else {
				H[0] = v8[0] - qe, max_H = H[0], max_t = 0; // special casing r==0
				need_t = ksw_extd2_need_max_t(ez, max_H, zdrop, zd_all);
			}
#ifdef KSW_H16
			if (H16) h_en0 = H16[en0], h_st0 = H16[st0];
			else
#endif
			h_en0 = H[en0], h_st0 = H[st0];
			// update ez
			if (en0 == tlen - 1 && h_en0 > ez->mte)
				ez->mte = h_en0, ez->mte_q = r - en0;
			if (r - st0 == qlen - 1 && h_st0 > ez->mqe)
				ez->mqe = h_st0, ez->mqe_t = st0;
			if (need_t && ksw_apply_zdrop(ez, 1, max_H, r, max_t, zdrop, e2)) break; // !need_t: it would return 0 and change nothing
			if (r == qlen + tlen - 2 && en0 == tlen - 1)
				ez->score = h_en0; // en0 == tlen - 1
			// early stopping in the extension-only mode; this block should not change the alignment
			if (flag & KSW_EZ_EXTZ_ONLY) {
				int32_t rH = last_max_H > max_H? last_max_H : max_H; // NB: if diagonal r is on the optimal path, r-1 or r+1 is not on the optimal path
				int32_t rq = qlen - (r - st0), rt = tlen - en0;
				// 3 conditions: hitting bottom of the DP matrix, hitting the right boundary of the matrix, and bracketing the diagonal
				int32_t rm = rq >= tlen - st0? tlen - st0 : rt >= qlen - (r - en0)? qlen - (r - en0) : tlen + qlen - 1 - r;
				if (rH + rm * max_sc + end_bonus < ez->max) break;
				// KSW_EZ_EXTZ_BOUND: every in-band cell's H is its diagonal predecessor's plus z, the core's z capped at mat[0] <= max_sc,
				// and that predecessor is in band too (same offset t-i), so a later cell, k diagonal moves after the cell of r or r-1
				// it leads back to, scores at most that cell's H + k * max_sc, k at most the cell's min(query, target bases left);
				// chains that reach row or column -1 after r-2 are bounded by the boundary's -min(q+l*e, q2+l*e2). When none of
				// that, plus end_bonus, reaches max, no later cell changes max, max_t, max_q, reach_end, or mqe/mqe_t if reach_end.
				// The bound is at least rH, so it is tried only below max; it is a scan of two anti-diagonals, so after one that does
				// not stop the loop, the next waits for it to fall by the excess at KSW_BOUND_DROP per anti-diagonal, unless max rises
				if (bound && rH + eb0 < (int32_t)ez->max) {
					if ((int32_t)ez->max > bmax) bmax = ez->max, bnext = r;
					if (r >= bnext) {
						int32_t lo1 = 0, hi1 = tlen - 1, B, x, g = q + (r + 1) * e < q2 + (r + 1) * e2? q + (r + 1) * e : q2 + (r + 1) * e2;
						if (lo1 < r - qlen) lo1 = r - qlen; // anti-diagonal r-1's band, as the loop's top computed it
						if (lo1 < (r-1-wr+1)>>1) lo1 = (r-1-wr+1)>>1;
						if (hi1 > r - 1) hi1 = r - 1;
						if (hi1 > (r-1+wl)>>1) hi1 = (r-1+wl)>>1;
#ifdef KSW_H16
						if (H16) B = ksw_extd2_bound16(H16, v8, st0, en0, lo1, hi1, qlen - 1 - r, tlen - 1, max_sc);
						else
#endif
						B = ksw_extd2_bound32(H, v8, st0, en0, lo1, hi1, qlen - 1 - r, tlen - 1, max_sc);
						if (r + 1 <= w && r + 1 <= tlen - 1) // later cells at offset t-i > r lead back to row -1
							x = -g + max_sc * (qlen < tlen - 1 - r? qlen : tlen - 1 - r), B = B > x? B : x;
						if (r + 1 <= w - 1 && r + 1 <= qlen - 1) // later cells at offset t-i < -r lead back to column -1
							x = -g + max_sc * (tlen < qlen - 1 - r? tlen : qlen - 1 - r), B = B > x? B : x;
						if (B + eb0 < (int32_t)ez->max) break;
						x = (B + eb0 - (int32_t)ez->max) / KSW_BOUND_DROP;
						bnext = r + (x < 1? 1 : x < KSW_BOUND_GAP? x : KSW_BOUND_GAP);
					}
				}
			}
			last_max_H = max_H; // need this for calculating rH
		} else { // find approximate max; Z-drop might be inaccurate, too.
			if (r > 0) {
				if (last_H0_t >= st0 && last_H0_t <= en0 && last_H0_t + 1 >= st0 && last_H0_t + 1 <= en0) {
					int32_t d0 = v8[last_H0_t];
					int32_t d1 = u8[last_H0_t + 1];
					if (d0 > d1) H0 += d0;
					else H0 += d1, ++last_H0_t;
				} else if (last_H0_t >= st0 && last_H0_t <= en0) {
					H0 += v8[last_H0_t];
				} else {
					++last_H0_t, H0 += u8[last_H0_t];
				}
			} else H0 = v8[0] - qe, last_H0_t = 0;
			if ((flag & KSW_EZ_APPROX_DROP) && ksw_apply_zdrop(ez, 1, H0, r, last_H0_t, zdrop, e2)) break;
			if (r == qlen + tlen - 2 && en0 == tlen - 1)
				ez->score = H0;
		}
		last_st = st, last_en = en;
		//for (t = st0; t <= en0; ++t) printf("(%d,%d)\t(%d,%d,%d,%d)\t%d\n", r, t, ((int8_t*)u)[t], ((int8_t*)v)[t], ((int8_t*)x)[t], ((int8_t*)y)[t], H[t]); // for debugging
	}
	if (!ksw_scr_has(sc, KSW_SCR_MEM, mem)) kfree(km, mem);
	if (!approx_max && !ksw_scr_has(sc, KSW_SCR_H, H)) kfree(km, H);
#ifdef KSW_H16
	if (!ksw_scr_has(sc, KSW_SCR_H16, H16)) kfree(km, H16);
#endif
	if (with_cigar) { // backtrack
		int rev_cigar = !!(flag & KSW_EZ_REV_CIGAR);
		if (!ez->zdropped && !(flag&KSW_EZ_EXTZ_ONLY)) {
			ksw_backtrack(km, 1, rev_cigar, 0, (uint8_t*)p, off, off_end, n_col_*16, tlen-1, qlen-1, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
		} else if ((flag&KSW_EZ_EXTZ_ONLY) && ez->mqe + end_bonus > (int)ez->max) {
			ez->reach_end = 1;
			ksw_backtrack(km, 1, rev_cigar, 0, (uint8_t*)p, off, off_end, n_col_*16, ez->mqe_t, qlen-1, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
		} else if (ez->max_t >= 0 && ez->max_q >= 0) {
			ksw_backtrack(km, 1, rev_cigar, 0, (uint8_t*)p, off, off_end, n_col_*16, ez->max_t, ez->max_q, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
		}
		if (!ksw_scr_has(sc, KSW_SCR_MEM2, mem2)) kfree(km, mem2);
		if (!ksw_scr_has(sc, KSW_SCR_OFF, off)) kfree(km, off);
	}
#undef __dp_code_block1
#undef __dp_code_block2
#if KSW_EXTD2_BODY == 256
#undef __dp_code_block1_256
#undef __dp_code_block2_256
#endif
}
#endif // KSW_EXTD2_BODY
