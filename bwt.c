#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include "kommon.h"
#include "kalloc.h"
#include "bwt.h"

/********************
 * Basic operations *
 ********************/

static void bwt_gen_cnt_table(uint32_t cnt[256])
{
	int i, j;
	for (i = 0; i != 256; ++i) {
		uint32_t x = 0;
		for (j = 0; j != 4; ++j)
			x |= (((i&3) == j) + ((i>>2&3) == j) + ((i>>4&3) == j) + (i>>6 == j)) << (j<<3);
		cnt[i] = x;
	}
}

mb_bwt_t *mb_bwt_init(void)
{
	mb_bwt_t *bwt;
	bwt = kom_calloc(mb_bwt_t, 1);
	bwt->sa_bit = (uint32_t)-1;
	bwt_gen_cnt_table(bwt->cnt_table);
	return bwt;
}

void mb_bwt_destroy(mb_bwt_t *bwt)
{
	if (bwt == 0) return;
	free(bwt->pre); // pre is always heap-allocated (by mb_bwt_cache)
	free(bwt->pre2); // and pre2 by mb_bwt_cache2
	if (bwt->mmap) kom_munmap(bwt->mmap, bwt->mmap_len); // data/sa point into the mapped file
	else { free(bwt->sa); free(bwt->data); }
	free(bwt);
}

/******************
 * Encode raw BWT *
 ******************/

#define raw_B00(b, k) ((b)[(k)>>4]>>((~(k)&0xf)<<1)&3)

#define BWT_CNT_SHIFT 56
#define BWT_CNT_MASK ((1ULL<<BWT_CNT_SHIFT) - 1)

static uint64_t mb_bwt_data_len(uint64_t len)
{
	uint64_t bwt_len, occ_len;
	bwt_len = (len + 127) / 128 * 4;
	occ_len = ((len + 127) / 128 + 1) * 4; // +1 for the final counts
	return bwt_len + occ_len;
}

/* BWT layout. Each block consists of u64[4]+u32[8], 64 bytes in total. The
 * lower 56 bits of each u64[4] (see BWT_CNT_SHIFT) store the accumulative
 * count of A/C/G/T bases. The higher 8 bits store the count of A/C/G/T in the
 * next 64nt. u32[8] keeps a BWT substring of 128nt in length. Because it
 * follows little endian, it can also be considered as u64[4] etc.
 */
mb_bwt_t *mb_bwt_init_from_raw(int is_byte, const void *raw_, uint64_t len, uint64_t primary)
{
	uint64_t c[4], x[4], i, k, *last_c = 0;
	mb_bwt_t *bwt;
	const uint32_t *raw32 = is_byte? 0 : (const uint32_t*)raw_;
	const uint8_t *raw8 = is_byte? (const uint8_t*)raw_ : 0;

	bwt = mb_bwt_init();
	bwt->primary = primary;
	bwt->seq_len = len;
	bwt->data_len = mb_bwt_data_len(len);
	bwt->data = kom_calloc(uint64_t, bwt->data_len);
	if (len == 0) return bwt; // nothing to encode; avoid the last-block overflow

	memset(c, 0, 32);
	for (i = k = 0; i < len; ++i) {
		uint8_t j, a = is_byte? raw8[i]&3 : raw_B00(raw32, i);
		if ((i & 0x7f) == 0) { // hard coded: each block encodes 128bp
			if (i > 0) {
				memcpy(&bwt->data[k], x, 32);
				k += 4;
			}
			last_c = &bwt->data[k];
			memcpy(&bwt->data[k], c, 32);
			k += 4;
			memset(x, 0, 32);
		} else if ((i & 0x3f) == 0 && last_c) {
			for (j = 0; j < 4; ++j)
				last_c[j] |= (c[j] - last_c[j]) << BWT_CNT_SHIFT;
		}
		++c[a];
		x[(i&0x7f)>>5] |= (uint64_t)a << ((i&0x1f)<<1); // little endian
	}
	// the last block
	memcpy(&bwt->data[k], x, 32);
	k += 4;
	memcpy(&bwt->data[k], c, 32);
	k += 4;
	assert(k == bwt->data_len);
	for (i = 0, bwt->L2[0] = 0; i < 4; ++i)
		bwt->L2[i+1] = bwt->L2[i] + c[i];
	assert(bwt->L2[4] == len);
	return bwt;
}

/********
 * Rank *
 ********/

#define bwt_block(b, k) ((b)->data + ((k)>>7<<3))

// retrieve a character from the $-removed BWT string. Note that mb_bwt_t::data is
// not exactly the BWT string and therefore this macro is called bwt_B0 instead of bwt_B
#define bwt_B0(b, k) ((b)->data[((k)>>7<<3) + 4 + (((k)&127)>>5)] >> (((k)&31)<<1) & 3)

/* always_inline: gcc -O3 marks a helper whose only effect is a prefetch as
 * pure and deletes every call to it before inlining, prefetch included */
static inline __attribute__((always_inline)) void mb_bwt_block_prefetch(const mb_bwt_t *bwt, uint64_t k)
{
	if (k > 0) __builtin_prefetch(bwt_block(bwt, k - 1 - (k - 1 >= bwt->primary)));
}

static inline int rank_aux1(uint64_t y, uint8_t c)
{
	// reduce nucleotide counting to bits counting
	y = ((c&2)? y : ~y) >> 1 & ((c&1)? y : ~y) & 0x5555555555555555ull;
	// count the number of 1s in y
	#if 0
	y = (y & 0x3333333333333333ull) + (y >> 2 & 0x3333333333333333ull);
	return ((y + (y >> 4)) & 0xf0f0f0f0f0f0f0full) * 0x101010101010101ull >> 56;
	#else
	return __builtin_popcountll(y);
	#endif
}

uint64_t mb_bwt_rank11(const mb_bwt_t *bwt, uint64_t k, uint8_t c)
{
	const uint64_t *p, *end;
	uint64_t n, mask;
	if (k == 0) return 0;
	if (k == bwt->seq_len + 1) return bwt->L2[c+1] - bwt->L2[c];
	--k;
	k -= (k >= bwt->primary); // because $ is not in bwt
	mask = (k&0x7f) >= 64? (1ULL << (64 - BWT_CNT_SHIFT)) - 1 : 0;
	p = bwt_block(bwt, k);
	n = (p[c] & BWT_CNT_MASK) + ((p[c] >> BWT_CNT_SHIFT) & mask);
	p += 4; // p points to 2-bit encoded BWT
	end = p + ((k&0x7f) >> 5);
	p += (k&0x7f) >= 64? 2 : 0;
//	for (; p < end; ++p) n += rank_aux1(*p, c); // we go through this loop 0 or 1 time
	if (p < end) n += rank_aux1(*p, c), ++p;
	n += rank_aux1(*p << ((~k&0x1f) << 1), c);
	if (c == 0) n -= ~k&0x1f; // "A" may be overcounted due to the shift above; this line corrects that
	return n;
}

#define BWT_EVEN 0x5555555555555555ULL

/* Counts of A/C/G/T in the $-removed BWT string up to and including position k:
 * the block's cumulative counts, plus the 8-bit counts of the block's first 64nt
 * when k is in the second half, plus the counts over the 1-64nt of k's half up
 * to k. The latter come from masked popcounts over the half's two u64 words
 * (little endian; position i in bits 2i and 2i+1): the word holding k is shifted
 * left so that positions after k drop out and zeros (A) enter at the bottom, so
 * C, G and T need no mask and A is the rest of the k&0x3f+1 positions. Both words
 * are folded into one u64 per popcount: word 0 on the even bits, word 1 on the odd. */
static inline void rank_half(const mb_bwt_t *bwt, uint64_t k, uint64_t cnt[4])
{
	const uint64_t *p = bwt_block(bwt, k), *q = p + 4 + (k>>5&2); // q: the 64nt half holding k
	uint64_t mask = (k&0x7f) >= 64? (1ULL << (64 - BWT_CNT_SHIFT)) - 1 : 0;
	uint64_t w = -(k>>5&1); // all ones when k is in the half's second word; no branch
	uint64_t y0, y1, ct, gt, t;
	uint32_t sh = (~k&0x1f) << 1, n1, n2, n3;
	y0 = (q[0] & w) | (q[0] << sh & ~w);
	y1 = q[1] << sh & w;
	ct = (y0 & BWT_EVEN) | (y1 & BWT_EVEN) << 1; // C or T: low bit set
	gt = (y0 >> 1 & BWT_EVEN) | (y1 & ~BWT_EVEN); // G or T: high bit set
	t = (y0 & y0 >> 1 & BWT_EVEN) | (y1 & y1 >> 1 & BWT_EVEN) << 1; // T: both bits set
	n3 = __builtin_popcountll(t);
	n1 = __builtin_popcountll(ct) - n3;
	n2 = __builtin_popcountll(gt) - n3;
	cnt[0] = (p[0] & BWT_CNT_MASK) + ((p[0] >> BWT_CNT_SHIFT) & mask) + ((k&0x3f) + 1 - n1 - n2 - n3);
	cnt[1] = (p[1] & BWT_CNT_MASK) + ((p[1] >> BWT_CNT_SHIFT) & mask) + n1;
	cnt[2] = (p[2] & BWT_CNT_MASK) + ((p[2] >> BWT_CNT_SHIFT) & mask) + n2;
	cnt[3] = (p[3] & BWT_CNT_MASK) + ((p[3] >> BWT_CNT_SHIFT) & mask) + n3;
}

static inline uint32_t rank_aux4(const mb_bwt_t *bwt, uint32_t x)
{
	return bwt->cnt_table[x&0xff] + bwt->cnt_table[x>>8&0xff] + bwt->cnt_table[x>>16&0xff] + bwt->cnt_table[x>>24];
}

/* cntl of mb_bwt_rank2a() when k and l are in one block, counted from k1:
 * the block's counts at k1, plus l's symbols counted word by word from the
 * start of k1's 64-row half. It gives rank_half(l1) whenever k1 <= l1 < seq_len,
 * so it is only called otherwise: on l < k, an interval of negative size, it
 * reads l1's offset in k1's word when l1's word comes before k1's (the loops
 * run zero times), and past the text (rows the BWT does not hold) it counts the
 * first half's words where rank_half() reads the block's stored counts. Such
 * intervals arise from the forward extension on a single-strand index (genbwt
 * -1), whose x[1] is not the reverse complement's interval, and the seeds of
 * such an index depend on how this function counts them. */
static __attribute__((noinline, cold)) void rank_l_from_k_half(const mb_bwt_t *bwt, uint64_t k1, uint64_t l1, uint64_t cntl[4])
{
	const uint64_t *p = bwt_block(bwt, k1);
	uint64_t mask = (k1&0x7f) >= 64? (1ULL << (64 - BWT_CNT_SHIFT)) - 1 : 0;
	const uint32_t *q = (const uint32_t*)(p + 4), *endk = q + ((k1&0x7f) >> 4), *endl = q + ((l1&0x7f) >> 4);
	uint32_t y = 0, tmp;
	q += (k1&0x7f) >= 64? 4 : 0;
	for (; q < endk; ++q) y += rank_aux4(bwt, *q); // NB: this assumes little endian
	for (; q < endl; ++q) y += rank_aux4(bwt, *q);
	tmp = *q << ((~l1&0xf) << 1);
	y += rank_aux4(bwt, tmp) - (~l1&0xf);
	cntl[0] = (p[0] & BWT_CNT_MASK) + ((p[0] >> BWT_CNT_SHIFT) & mask) + (y&0xff);
	cntl[1] = (p[1] & BWT_CNT_MASK) + ((p[1] >> BWT_CNT_SHIFT) & mask) + (y>>8&0xff);
	cntl[2] = (p[2] & BWT_CNT_MASK) + ((p[2] >> BWT_CNT_SHIFT) & mask) + (y>>16&0xff);
	cntl[3] = (p[3] & BWT_CNT_MASK) + ((p[3] >> BWT_CNT_SHIFT) & mask) + (y>>24);
}

void mb_bwt_rank1a(const mb_bwt_t *bwt, uint64_t k, uint64_t cnt[4])
{
	if (k == 0) {
		memset(cnt, 0, 4 * sizeof(uint64_t));
		return;
	}
	--k;
	k -= (k >= bwt->primary); // because $ is not in bwt
	rank_half(bwt, k, cnt);
}

void mb_bwt_rank2a(const mb_bwt_t *bwt, uint64_t k, uint64_t l, uint64_t cntk[4], uint64_t cntl[4])
{
	uint64_t k1 = k - 1, l1 = l - 1;
	k1 -= (k1 >= bwt->primary);
	l1 -= (l1 >= bwt->primary);
	mb_bwt_block_prefetch(bwt, k);
	if (k1>>7 != l1>>7 || k == 0 || l == 0) {
		mb_bwt_block_prefetch(bwt, l);
		mb_bwt_rank1a(bwt, k, cntk);
		mb_bwt_rank1a(bwt, l, cntl);
	} else if (l - k == 1) { // we can use a simpler procedure
		uint64_t z = k - (k > bwt->primary);
		mb_bwt_rank1a(bwt, k, cntk);
		memcpy(cntl, cntk, 4 * sizeof(uint64_t));
		++cntl[bwt_B0(bwt, z)];
	} else { // k and l in the same block
		rank_half(bwt, k1, cntk);
		if (l1 < k1 || l1 >= bwt->seq_len) rank_l_from_k_half(bwt, k1, l1, cntl);
		else rank_half(bwt, l1, cntl);
	}
}

/*********************
 * Bidirectional BWT *
 *********************/

void mb_bwt_extend(const mb_bwt_t *bwt, const mb_sai_t *ik, mb_sai_t ok[4], int is_back)
{
	uint64_t tk[4], tl[4];
	int i;
	mb_bwt_rank2a(bwt, ik->x[!is_back], ik->x[!is_back] + ik->size, tk, tl);
	for (i = 0; i != 4; ++i) {
		ok[i].x[!is_back] = bwt->L2[i] + 1 + tk[i]; // +1 for the missing sentinel
		ok[i].size = (tl[i] -= tk[i]);
	}
	ok[3].x[is_back] = ik->x[is_back] + (ik->x[!is_back] <= bwt->primary && ik->x[!is_back] + ik->size > bwt->primary);
	ok[2].x[is_back] = ok[3].x[is_back] + tl[3];
	ok[1].x[is_back] = ok[2].x[is_back] + tl[2];
	ok[0].x[is_back] = ok[1].x[is_back] + tl[1];
}

// backward search from pos
static int64_t mb_bwt_back(const mb_bwt_t *f, const uint8_t *q, int64_t st, int64_t pos, int64_t min_occ, mb_sai_t *p)
{
	int64_t i = pos - 1;
	mb_sai_t ok[4];
	assert(q[pos] < 4); // the backward pass never involves N
	if (f->pre2 && pos - st >= f->pre2_len) { // the longer k-mer; if it jumps too far, the pre_len-mer ending at pos
		uint64_t z = 0, l = 0;
		for (i = pos; l < f->pre2_len; --i, ++l)
			z = z << 2 | q[i];
		*p = f->pre2[z];
		if (p->size < min_occ) {
			i += f->pre2_len - f->pre_len;
			*p = f->pre[z >> 2 * (f->pre2_len - f->pre_len)];
		}
	} else if (f->pre && pos - st >= f->pre_len) { // then use precomputed k-mer index instead of base-by-base extension
		uint64_t z = 0, l = 0;
		for (i = pos; l < f->pre_len; --i, ++l) // get the k-mer
			z = z << 2 | q[i];                  // NB: this loop doesn't check N
		assert(z < 1<<f->pre_len*2);
		*p = f->pre[z];
	} else p->size = 0;
	if (p->size < min_occ) { // then we need to use the standard procedure
		mb_bwt_set_intv(f, q[pos], p);
		i = pos - 1;
	}
	for (; i >= st; --i) { // backward extension
		int c = q[i];
		if (c > 3) break;
		mb_bwt_extend(f, p, ok, 1);
		if (ok[c].size < min_occ) break;
		*p = ok[c];
	}
	return i;
}

// find super MEMs (SMEMs). See ropebwt3
int64_t mb_bwt_smem(const mb_bwt_t *f, uint32_t len, const uint8_t *q, int64_t x, int64_t min_len, int64_t min_occ, mb_sai_t *p)
{
	int64_t i, j, xn;
	mb_sai_t ik, ok[4];

	assert(len <= INT32_MAX); // this can be relaxed if we define a new struct for mem
	p->size = ik.size = 0;
	if (len - x < min_len) return len;
	for (i = x, xn = -1; i < x + min_len; ++i) // find the last N in [x,x+min_len)
		if (q[i] > 3) xn = i;
	if (xn >= 0) return xn + 1;
	i = mb_bwt_back(f, q, x, x + min_len - 1, min_occ, &ik);
	if (i >= x) return i + 1; // no MEM found
	for (j = x + min_len; j < len; ++j) { // forward extension
		int c = 3 - q[j];
		if (q[j] > 3) break;
		mb_bwt_extend(f, &ik, ok, 0);
		if (ok[c].size < min_occ) break;
		ik = ok[c];
	}
	*p = ik;
	p->info = (uint64_t)x<<32 | j;
	if (j == len) return len; // reaching end; no need to do another round
	i = q[j] > 3? j : mb_bwt_back(f, q, x + 1, j, min_occ, &ik);
	return i + 1;
}

/**************
 * Batch SMEM *
 **************/

static inline uint64_t bwt_invPsi(const mb_bwt_t *bwt, uint64_t k) // compute inverse CSA
{
	uint64_t x = k - (k > bwt->primary);
	int c = bwt_B0(bwt, x);
	x = bwt->L2[c] + 1 + mb_bwt_rank11(bwt, k, c); // +1 to account for the sentinel
	return k == bwt->primary? 0 : x;
}

typedef struct { // a simplified version of kdq
	int32_t front, count, cap;
	int32_t *a;
} tiny_queue_t;

static void tq_init(void *km, tiny_queue_t *q, int32_t n)
{
	q->cap = n;
	kom_roundup32(q->cap);
	q->a = Kcalloc(km, int32_t, q->cap);
	q->front = q->count = 0;
}

static inline void tq_push(tiny_queue_t *q, int32_t x)
{
	q->a[((q->count++) + q->front) & (q->cap - 1)] = x;
}

static inline int32_t tq_shift(tiny_queue_t *q)
{
	int32_t x;
	if (q->count == 0) return -1;
	x = q->a[q->front++];
	q->front &= q->cap - 1;
	--q->count;
	return x;
}

/* rank_half() for se_extend1(), with the four counts in scalars. A private
 * always_inline copy: another caller of the shared rank code changes how gcc
 * inlines it into its other callers */
static inline __attribute__((always_inline)) void se_rank_half(const mb_bwt_t *bwt, uint64_t k, uint64_t *c0, uint64_t *c1, uint64_t *c2, uint64_t *c3)
{
	const uint64_t *p = bwt_block(bwt, k), *q = p + 4 + (k>>5&2);
	uint64_t mask = (k&0x7f) >= 64? (1ULL << (64 - BWT_CNT_SHIFT)) - 1 : 0;
	uint64_t w = -(k>>5&1);
	uint64_t y0, y1, ct, gt, t;
	uint32_t sh = (~k&0x1f) << 1, n1, n2, n3;
	y0 = (q[0] & w) | (q[0] << sh & ~w);
	y1 = q[1] << sh & w;
	ct = (y0 & BWT_EVEN) | (y1 & BWT_EVEN) << 1;
	gt = (y0 >> 1 & BWT_EVEN) | (y1 & ~BWT_EVEN);
	t = (y0 & y0 >> 1 & BWT_EVEN) | (y1 & y1 >> 1 & BWT_EVEN) << 1;
	n3 = __builtin_popcountll(t);
	n1 = __builtin_popcountll(ct) - n3;
	n2 = __builtin_popcountll(gt) - n3;
	*c0 = (p[0] & BWT_CNT_MASK) + ((p[0] >> BWT_CNT_SHIFT) & mask) + ((k&0x3f) + 1 - n1 - n2 - n3);
	*c1 = (p[1] & BWT_CNT_MASK) + ((p[1] >> BWT_CNT_SHIFT) & mask) + n1;
	*c2 = (p[2] & BWT_CNT_MASK) + ((p[2] >> BWT_CNT_SHIFT) & mask) + n2;
	*c3 = (p[3] & BWT_CNT_MASK) + ((p[3] >> BWT_CNT_SHIFT) & mask) + n3;
}

/* rank_half(bwt, k, cnt) reduced to cnt[c] (*n) and the sum of cnt[] above c
 * (*u), selected by masks: c is the next base, so a branch on it mispredicts */
static inline __attribute__((always_inline)) void se_rank_c(const mb_bwt_t *bwt, uint64_t k, int32_t c, uint64_t *n, uint64_t *u)
{
	uint64_t c0, c1, c2, c3;
	uint64_t e0 = -(uint64_t)(c == 0), e1 = -(uint64_t)(c == 1), e2 = -(uint64_t)(c == 2), e3 = -(uint64_t)(c == 3);
	se_rank_half(bwt, k, &c0, &c1, &c2, &c3);
	*n = (c0 & e0) | (c1 & e1) | (c2 & e2) | (c3 & e3);
	*u = (c3 & ~e3) + (c2 & (e0 | e1)) + (c1 & e0);
}

/* The queue's extension of p by the symbol c: mb_bwt_extend(bwt, p, ok, is_back)
 * followed by `if (ok[c].size >= min_occ) *p = ok[c]`, forming ok[c] only, with
 * mb_bwt_rank2a()'s counts in registers, so that no load reads back a store of
 * another width. ok[c].size is cl - ck and the sizes above c sum to ul - uk,
 * both modulo 2^64 as in mb_bwt_extend(). Returns whether p was extended;
 * p->info is left as it is (mb_bwt_extend() never writes ok[].info, and the
 * queue sets p->info before reading it) */
static inline __attribute__((always_inline)) int se_extend1(const mb_bwt_t *bwt, mb_sai_t *p, int32_t c, int is_back, int32_t min_occ)
{
	uint64_t k = p->x[!is_back], l = k + p->size, k1 = k - 1, l1 = l - 1;
	uint64_t ck, uk, cl, ul, sz;
	int same; // as mb_bwt_rank2a(): k and l in one block, neither 0
	k1 -= (k1 >= bwt->primary);
	l1 -= (l1 >= bwt->primary);
	same = !(k1>>7 != l1>>7 || k == 0 || l == 0);
	mb_bwt_block_prefetch(bwt, k);
	if (!same) mb_bwt_block_prefetch(bwt, l);
	if (k == 0) ck = uk = 0;
	else se_rank_c(bwt, k1, c, &ck, &uk);
	if (same && l - k == 1) { // the counts at k plus the symbol at row k
		int32_t b = bwt_B0(bwt, k - (k > bwt->primary));
		cl = ck + (b == c), ul = uk + (b > c);
	} else if (l == 0) cl = ul = 0;
	else if (same && (l1 < k1 || l1 >= bwt->seq_len)) { // as mb_bwt_rank2a(); never on an interval within the text
		uint64_t t[4];
		int32_t i;
		rank_l_from_k_half(bwt, k1, l1, t);
		cl = t[c];
		for (i = c + 1, ul = 0; i < 4; ++i) ul += t[i];
	} else se_rank_c(bwt, l1, c, &cl, &ul);
	sz = cl - ck;
	if (sz < min_occ) return 0;
	p->x[!is_back] = bwt->L2[c] + 1 + ck; // +1 for the missing sentinel
	p->x[is_back] += (k <= bwt->primary && l > bwt->primary) + (ul - uk); // plus the sizes of the symbols above c
	p->size = sz;
	return 1;
}

static inline __attribute__((always_inline)) void se_one_step_back(const mb_bwt_t *bwt, mb_smem_entry_t *s)
{
	int32_t c = s->q[s->i];
	assert(c < 4); // shouldn't happen
	if (!se_extend1(bwt, &s->p, c, 1, s->min_occ)) { // move back to stage1
		s->x = s->i + 1;
		s->stage = 1;
	} else { // stay in the two backward stages
		s->i--;
		mb_bwt_block_prefetch(bwt, s->p.x[0]); // prefetch for the next backward iteration
		mb_bwt_block_prefetch(bwt, s->p.x[0] + s->p.size);
	}
}

/* T[t] of the BWT text, the 2-bit forward strand in pac followed by its
 * reverse complement (seq_len = 2*tot_len), compared with q[i..en) from
 * T[t]; returns the first i at which they differ (q[i] is N or T has ended)
 * or en */
static inline int32_t se_text_match(const uint64_t *pac, uint64_t seq_len, const uint8_t *q, int32_t i, int32_t en, uint64_t t)
{
	uint64_t half = seq_len >> 1;
	for (; i < en && t < half; ++i, ++t)
		if (q[i] != (pac[t>>5] >> ((t&31)<<1) & 3)) return i;
	for (; i < en && t < seq_len; ++i, ++t) {
		uint64_t u = seq_len - 1 - t;
		if (q[i] != 3 - (pac[u>>5] >> ((u&31)<<1) & 3)) return i;
	}
	return i;
}

static inline void se_text_prefetch(const uint64_t *pac, uint64_t seq_len, uint64_t t)
{
	uint64_t half = seq_len >> 1;
	if (t < half) __builtin_prefetch(&pac[t>>5]);
	else if (t < seq_len) __builtin_prefetch(&pac[(seq_len - 1 - t)>>5]);
}

/* Enter or stay in the forward pass at s->i with s->p. With the text in pac,
 * an interval of size 1 extends by a base iff its one occurrence is followed
 * by that base, so it goes to the SA walk (stages 7 and 8, ->p.x[1] holding
 * the row and ->p.info the steps) and then reads the text (stage 9). */
static inline void se_forward(const mb_bwt_t *bwt, const uint64_t *pac, mb_smem_entry_t *s)
{
	if (pac && s->p.size == 1 && s->min_occ == 1 && s->i < s->en && s->q[s->i] < 4) {
		uint64_t k = s->p.x[0];
		s->p.x[1] = k, s->p.info = 0;
		if ((k & ((1ULL<<bwt->sa_bit) - 1)) == 0) {
			__builtin_prefetch(&bwt->sa[k >> bwt->sa_bit]);
			s->stage = 8;
		} else {
			mb_bwt_block_prefetch(bwt, k);
			s->stage = 7;
		}
	} else {
		mb_bwt_block_prefetch(bwt, s->p.x[1]); // prefetch for the forward pass
		mb_bwt_block_prefetch(bwt, s->p.x[1] + s->p.size);
		s->stage = 4;
	}
}

#define MB_KMER2 0x80000000U // in mb_smem_entry_t::kmer: a key of bwt->pre2, not of bwt->pre

/* Reads the k-mer of s->q ending at s->i into s->kmer and prefetches its table
 * entry: the pre2_len-mer (flagged MB_KMER2) when it fits in the n bases
 * available, else the pre_len-mer, which the caller checked fits; s->i ends
 * before the k-mer. The backward pass never meets N here. */
static inline void se_kmer_get(const mb_bwt_t *bwt, mb_smem_entry_t *s, int32_t n)
{
	int32_t i;
	if (bwt->pre2 && n >= bwt->pre2_len) {
		for (i = 0, s->kmer = 0; i < bwt->pre2_len; ++i, s->i--)
			s->kmer = s->kmer << 2 | s->q[s->i];
		__builtin_prefetch(&bwt->pre2[s->kmer]);
		s->kmer |= MB_KMER2;
	} else {
		for (i = 0, s->kmer = 0; i < bwt->pre_len; ++i, s->i--)
			s->kmer = s->kmer << 2 | s->q[s->i];
		__builtin_prefetch(&bwt->pre[s->kmer]);
	}
}

// the forward pass can't extend s->p at s->i < s->en: save the SMEM and start the second backward pass
static inline void se_forward_end(void *km, const mb_bwt_t *bwt, mb_smem_entry_t *s)
{
	s->p.info = (uint64_t)s->x << 32 | s->i;
	Kgrow(km, mb_sai_t, s->v->a, s->v->n, s->v->m);
	s->v->a[s->v->n++] = s->p; // save the interval
	if (s->q[s->i] > 3) { // if N, skip it and move back to stage 1
		s->x = s->i + 1;
		s->stage = 1;
	} else if (bwt->pre && s->i - s->x - 1 >= bwt->pre_len) { // get k-mer
		se_kmer_get(bwt, s, s->i - s->x - 1);
		s->stage = 5;
	} else { // skip stage 5
		mb_bwt_set_intv(bwt, s->q[s->i--], &s->p);
		s->stage = 6;
	}
}

static inline void smem_entry_init(void *km, mb_smem_entry_t *s)
{
	s->stage = 1;
	s->x = s->st;
	if (s->v->m < 64) { // preallocate to avoid frequent krealloc(), which can be slow
		s->v->m = 64;
		s->v->a = Krealloc(km, mb_sai_t, s->v->a, s->v->m);
	}
}

/* One step of the entry s of mb_bwt_smem_feed(); returns 0 when s has finished.
 * always_inline: the step is the batch SMEM's inner loop */
static inline __attribute__((always_inline)) int smem_step(void *km, const mb_bwt_t *bwt, const uint64_t *pac, mb_smem_entry_t *s)
{
	if (s->stage == 1) { // set interval for the first backward pass in smem; require ->x
		int32_t i, xn;
		if (s->en - s->x < s->min_len)
			return 0; // the entry has finished
		for (i = s->x, xn = -1; i < s->x + s->min_len; ++i) // find the position of the last N
			if (s->q[i] > 3) xn = i;
		if (xn >= 0) { // skip N and stay in stage 1
			s->x = xn + 1;
		} else {
			s->i = s->x + s->min_len - 1;
			if (bwt->pre && s->min_len >= bwt->pre_len) { // get k-mer for prefetch
				se_kmer_get(bwt, s, s->min_len); // backward pass shouldn't meet N
				s->stage = 2;
			} else { // skip stage 2
				mb_bwt_set_intv(bwt, s->q[s->i--], &s->p);
				s->stage = 3;
			}
		}
	} else if (s->stage == 2 || s->stage == 5) { // k-mer lookup
		if (s->kmer & MB_KMER2) { // the longer k-mer; if it jumps too far, the pre_len-mer ending at the same base
			uint32_t z = s->kmer & ~MB_KMER2;
			s->p = bwt->pre2[z];
			if (s->p.size < s->min_occ) {
				s->i += bwt->pre2_len - bwt->pre_len;
				s->p = bwt->pre[z >> 2 * (bwt->pre2_len - bwt->pre_len)];
			}
		} else s->p = bwt->pre[s->kmer];
		if (s->p.size < s->min_occ) { // jumped too far with the k-mer cache; revert
			s->i += bwt->pre_len;
			mb_bwt_set_intv(bwt, s->q[s->i--], &s->p);
		}
		mb_bwt_block_prefetch(bwt, s->p.x[0]);
		mb_bwt_block_prefetch(bwt, s->p.x[0] + s->p.size);
		s->stage++;
	} else if (s->stage == 3 || s->stage == 6) { // first (3) or second (6) backward pass; require ->{i,p}
		if (s->i >= s->x + (s->stage == 6)) se_one_step_back(bwt, s);
		else if (s->stage == 3) { // move to the next stage
			s->i = s->x + s->min_len;
			se_forward(bwt, pac, s);
		} else {
			s->x = s->i + 1;
			s->stage = 1;
		}
	} else if (s->stage == 4) { // forward pass; require ->{i,p}
		if (s->i == s->en) {
			s->p.info = (uint64_t)s->x << 32 | s->i;
			Kgrow(km, mb_sai_t, s->v->a, s->v->n, s->v->m);
			s->v->a[s->v->n++] = s->p; // save the interval
			return 0; // the entry has finished
		} else {
			int32_t c = 3 - (int32_t)s->q[s->i];
			if (c >= 0 && se_extend1(bwt, &s->p, c, 0, s->min_occ)) { // stay in stage 4
				s->i++;
				se_forward(bwt, pac, s);
			} else se_forward_end(km, bwt, s);
		}
	} else if (s->stage == 7) { // SA walk of the size-1 interval's row, one step per turn
		uint64_t k = bwt_invPsi(bwt, s->p.x[1]);
		s->p.x[1] = k, ++s->p.info;
		if ((k & ((1ULL<<bwt->sa_bit) - 1)) == 0) {
			__builtin_prefetch(&bwt->sa[k >> bwt->sa_bit]);
			s->stage = 8;
		} else mb_bwt_block_prefetch(bwt, k);
	} else if (s->stage == 8) { // sampled SA value: the text position of x[0], as mb_bwt_sa_batch() gets it
		uint64_t pos = s->p.info + bwt->sa[s->p.x[1] >> bwt->sa_bit];
		s->p.x[1] = pos | MB_SA_KNOWN;
		se_text_prefetch(pac, bwt->seq_len, pos + (s->i - s->x));
		s->stage = 9;
	} else if (s->stage == 9) { // forward pass on the text
		uint64_t pos = s->p.x[1] & ~MB_SA_KNOWN;
		s->i = se_text_match(pac, bwt->seq_len, s->q, s->i, s->en, pos + (s->i - s->x));
		if (s->i == s->en) {
			s->p.info = (uint64_t)s->x << 32 | s->i;
			Kgrow(km, mb_sai_t, s->v->a, s->v->n, s->v->m);
			s->v->a[s->v->n++] = s->p; // save the interval
			return 0; // the entry has finished
		}
		se_forward_end(km, bwt, s);
	}
	return 1;
}

/* The batch SMEM queue over the n slots of a[]. Without feed, a[] holds n
 * entries and the call returns when all have finished. With feed, a slot is
 * filled by feed(data, slot, 0) at the start, and when its entry finishes by
 * feed(data, slot, 1), so that the queue stays full while jobs are pending;
 * with pac, as mb_bwt_smem_batch_pac(). */
void mb_bwt_smem_feed(void *km, const mb_bwt_t *bwt, const uint64_t *pac, int32_t n, mb_smem_entry_t *a, mb_smem_feed_f feed, void *data)
{
	int32_t i, n_free = 0, *free_slot = 0;
	tiny_queue_t tq;

	// initialize
	tq_init(km, &tq, n);
	if (feed) free_slot = Kmalloc(km, int32_t, n);
	for (i = 0; i < n; ++i) {
		if (feed && !feed(data, &a[i], 0)) {
			free_slot[n_free++] = i;
			continue;
		}
		tq_push(&tq, i);
		smem_entry_init(km, &a[i]);
	}

	// core loop
	while (tq.count > 0) {
		int32_t idx = tq_shift(&tq);
		if (smem_step(km, bwt, pac, &a[idx])) {
			tq_push(&tq, idx);
		} else if (feed) { // the entry has finished: refill its slot, and idle slots if feed() has more
			if (feed(data, &a[idx], 1)) {
				tq_push(&tq, idx);
				smem_entry_init(km, &a[idx]);
				while (n_free > 0 && feed(data, &a[free_slot[n_free - 1]], 0)) {
					idx = free_slot[--n_free];
					tq_push(&tq, idx);
					smem_entry_init(km, &a[idx]);
				}
			} else free_slot[n_free++] = idx;
		}
	}
	kfree(km, free_slot);
	kfree(km, tq.a);
}

void mb_bwt_smem_batch(void *km, const mb_bwt_t *bwt, int32_t n, mb_smem_entry_t *a)
{
	mb_bwt_smem_feed(km, bwt, 0, n, a, 0, 0);
}

/* With pac non-NULL, the 2-bit forward strand of the BWT text (as l2b_t::pac;
 * bwt->seq_len must be twice its length), a size-1 interval in the forward
 * pass with min_occ 1 is extended by reading the text at its occurrence, and
 * its x[1] is saved as that occurrence's text position | MB_SA_KNOWN instead
 * of the SA row of the reverse complement. x[0], size and info are unchanged. */
void mb_bwt_smem_batch_pac(void *km, const mb_bwt_t *bwt, const uint64_t *pac, int32_t n, mb_smem_entry_t *a)
{
	mb_bwt_smem_feed(km, bwt, pac, n, a, 0, 0);
}

/***************************
 * Suffix array operations *
 ***************************/

// bwt->bwt and bwt->occ must be precalculated
void mb_bwt_gen_sa(mb_bwt_t *bwt, uint32_t sa_bit)
{
	uint64_t isa, sa, i, mask; // S(isa) = sa

	assert(bwt->data);
	if (bwt->sa) free(bwt->sa);
	bwt->sa_bit = sa_bit;
	bwt->n_sa = (bwt->seq_len + (1<<sa_bit)) >> sa_bit;
	bwt->sa = kom_calloc(uint64_t, bwt->n_sa);

	// calculate SA value
	isa = 0, sa = bwt->seq_len, mask = (1ULL<<sa_bit) - 1;
	for (i = 0; i < bwt->seq_len; ++i) {
		if ((isa & mask) == 0) bwt->sa[isa >> bwt->sa_bit] = sa;
		--sa;
		isa = bwt_invPsi(bwt, isa);
	}
	if ((isa & mask) == 0) bwt->sa[isa >> bwt->sa_bit] = sa;
	bwt->sa[0] = (uint64_t)-1; // before this line, bwt->sa[0] = bwt->seq_len
}

uint64_t mb_bwt_sa(const mb_bwt_t *bwt, uint64_t k)
{
	uint64_t sa = 0, mask = (1ULL<<bwt->sa_bit) - 1;
	while (k & mask) {
		++sa;
		k = bwt_invPsi(bwt, k);
	}
	// without setting bwt->sa[0] = -1, the following line should be
	// changed to (sa + bwt->sa[k/bwt->sa_intv]) % (bwt->seq_len + 1)
	return sa + bwt->sa[k >> bwt->sa_bit];
}

void mb_bwt_sa_batch(void *km, const mb_bwt_t *bwt, int64_t n, uint64_t *x)
{
	uint64_t mask = (1ULL<<bwt->sa_bit) - 1;
	int64_t i, step = 0, r;
	kom128_t *z;
	if (n <= 0) return;
	z = Kmalloc(km, kom128_t, n);
	for (i = 0, r = 0; i < n; ++i) {
		if (x[i] & MB_SA_KNOWN) { // already a text position
			x[i] &= ~MB_SA_KNOWN;
			continue;
		}
		z[r].x = x[i], z[r].y = i;
		if ((z[r].x & mask) == 0)
			__builtin_prefetch(&bwt->sa[z[r].x >> bwt->sa_bit]);
		else
			mb_bwt_block_prefetch(bwt, z[r].x);
		++r;
	}
	for (step = 0; r > 0; ++step) {
		int64_t r0 = r;
		for (i = 0, r = 0; i < r0; ++i) {
			if ((z[i].x & mask) == 0)
				x[z[i].y] = step + bwt->sa[z[i].x >> bwt->sa_bit];
			else z[r++] = z[i];
		}
		for (i = 0; i < r; ++i) {
			z[i].x = bwt_invPsi(bwt, z[i].x);
			if ((z[i].x & mask) == 0)
				__builtin_prefetch(&bwt->sa[z[i].x >> bwt->sa_bit]);
			else
				mb_bwt_block_prefetch(bwt, z[i].x);
		}
	}
	kfree(km, z);
}

/******************
 * k-mer counting *
 ******************/

typedef struct {
	mb_sai_t p;
	int d, c;
} count_pair64_t;

void mb_bwt_count_kmer(const mb_bwt_t *bwt, int32_t depth, mb_sai_t *s) // adapted from kount in ropebwt3
{
	count_pair64_t *p, stack[64];
	int32_t i, a, s_top = 0;
	uint8_t str[16];
	assert(depth <= 15);
	for (a = 0; a < 4; ++a) {
		p = &stack[s_top++];
		mb_bwt_set_intv(bwt, a, &p->p);
		p->d = 1, p->c = a;
	}
	while (s_top > 0) {
		count_pair64_t top = stack[--s_top];
		mb_sai_t ok[4];
		if (top.d > 0) str[depth - top.d] = top.c;
		mb_bwt_extend(bwt, &top.p, ok, 1);
		for (a = 0; a < 4; ++a) {
			str[depth - top.d - 1] = a;
			if (top.d != depth - 1) {
				p = &stack[s_top++];
				p->p = ok[a];
				p->d = top.d + 1;
				p->c = a;
			} else { // reaching the length; store in s[]
				uint64_t x = 0;
				for (i = 0; i < depth; ++i)
					x |= (uint64_t)str[i] << i * 2;
				s[x] = ok[a];
			}
		}
	}
}

void mb_bwt_cache(mb_bwt_t *bwt, int32_t len)
{
	if (bwt->pre) free(bwt->pre);
	free(bwt->pre2); // derived from pre
	bwt->pre2 = 0, bwt->pre2_len = 0;
	bwt->pre_len = len;
	bwt->pre = kom_calloc(mb_sai_t, 1 << len*2);
	mb_bwt_count_kmer(bwt, len, bwt->pre);
}

/* The pre2_len-mer table, pre2_len = pre_len + d with d 1 or 2, built by d
 * levels of backward extension from the pre_len-mer entries, exactly as the
 * backward pass extends an interval it took from pre: an entry is the interval
 * the pass reaches from the pre_len-mer, whatever its size, and its info is 0.
 * Prepending base c to the pattern of key z gives key z<<2|c. A pfor item
 * extends CACHE2_BLK pre_len-mers and then their children, each level
 * prefetching the BWT blocks of the interval CACHE2_DIST ahead. The prefetch is
 * a macro: gcc -O3 deletes calls to a helper whose only effect is a prefetch
 * unless it is always_inline (mb_bwt_block_prefetch()). */
#define CACHE2_BLK 256
#define CACHE2_DIST 16

#define cache2_prefetch(bwt, p) do { \
		uint64_t k_ = (p)->x[0], l_ = k_ + (p)->size; \
		if (k_ > 0) __builtin_prefetch(bwt_block(bwt, k_ - 1 - (k_ - 1 >= (bwt)->primary))); \
		if (l_ > 0) __builtin_prefetch(bwt_block(bwt, l_ - 1 - (l_ - 1 >= (bwt)->primary))); \
	} while (0)

typedef struct {
	const mb_bwt_t *bwt;
	mb_sai_t *pre2;
	int32_t d;
} cache2_aux_t;

static void cache2_worker(void *data, long j, int tid)
{
	const cache2_aux_t *a = (const cache2_aux_t*)data;
	const mb_bwt_t *bwt = a->bwt;
	const mb_sai_t *pre = bwt->pre;
	uint64_t z0 = (uint64_t)j * CACHE2_BLK, n = 1ULL << 2 * bwt->pre_len, z1, z, r, nm;
	mb_sai_t mid[CACHE2_BLK * 4], ok[4], *out;
	int c;
	(void)tid;
	z1 = z0 + CACHE2_BLK < n? z0 + CACHE2_BLK : n;
	out = a->d == 1? &a->pre2[z0 << 2] : mid;
	for (z = z0; z < z1 && z < z0 + CACHE2_DIST; ++z)
		cache2_prefetch(bwt, &pre[z]);
	for (z = z0; z < z1; ++z) {
		if (z + CACHE2_DIST < z1) cache2_prefetch(bwt, &pre[z + CACHE2_DIST]);
		mb_bwt_extend(bwt, &pre[z], ok, 1);
		for (c = 0; c < 4; ++c) {
			mb_sai_t *o = &out[(z - z0) << 2 | c];
			*o = ok[c], o->info = 0;
		}
	}
	if (a->d == 1) return;
	nm = (z1 - z0) << 2;
	for (r = 0; r < nm && r < CACHE2_DIST; ++r)
		cache2_prefetch(bwt, &mid[r]);
	for (r = 0; r < nm; ++r) {
		if (r + CACHE2_DIST < nm) cache2_prefetch(bwt, &mid[r + CACHE2_DIST]);
		mb_bwt_extend(bwt, &mid[r], ok, 1);
		for (c = 0; c < 4; ++c) {
			mb_sai_t *o = &a->pre2[((z0 << 2) + r) << 2 | c];
			*o = ok[c], o->info = 0;
		}
	}
}

void mb_bwt_cache2(mb_bwt_t *bwt, int32_t len, int n_threads, kom_pfor_f pfor)
{
	cache2_aux_t a;
	long n_blk, j;
	free(bwt->pre2);
	bwt->pre2 = 0, bwt->pre2_len = 0;
	if (bwt->pre == 0 || len <= (int32_t)bwt->pre_len || len > (int32_t)bwt->pre_len + 2 || len > 15)
		return; // 15: mb_smem_entry_t::kmer keeps bit 31 for MB_KMER2
	a.bwt = bwt, a.d = len - bwt->pre_len;
	a.pre2 = kom_malloc(mb_sai_t, 1ULL << 2 * len);
	if (a.pre2 == 0) return; // seeding is the same without it
	n_blk = ((1L << 2 * bwt->pre_len) + CACHE2_BLK - 1) / CACHE2_BLK;
	if (n_threads > 1 && pfor) pfor(n_threads, cache2_worker, &a, n_blk);
	else for (j = 0; j < n_blk; ++j) cache2_worker(&a, j, 0);
	bwt->pre2 = a.pre2, bwt->pre2_len = len;
}

/*************************
 * Read/write BWT and SA *
 *************************/

static uint64_t read_huge(FILE *fp, uint64_t size, void *a)
{ // Mac/Darwin has a bug when reading data longer than 2GB. This function fixes this issue by reading data in small chunks
	const int bufsize = 0x1000000; // 16M block
	uint64_t offset = 0;
	while (size) {
		int x = bufsize < size? bufsize : size;
		if ((x = fread(a + offset, 1, x, fp)) == 0) break;
		size -= x; offset += x;
	}
	return offset;
}

mb_bwt_t *mb_bwt_load_raw(const char *fn)
{
	mb_bwt_t *bwt;
	uint32_t *raw;
	uint64_t L2[5], primary, raw_size;
	FILE *fp;

	fp = fopen(fn, "rb");
	fseek(fp, 0, SEEK_END);
	raw_size = (ftell(fp) - sizeof(uint64_t) * 5) >> 2;
	raw = kom_calloc(uint32_t, raw_size);
	fseek(fp, 0, SEEK_SET);
	fread(&primary, sizeof(uint64_t), 1, fp);
	fread(L2 + 1, sizeof(uint64_t), 4, fp);
	L2[0] = 0;
	read_huge(fp, raw_size<<2, raw);
	fclose(fp);
	bwt = mb_bwt_init_from_raw(0, raw, L2[4], primary);
	free(raw);
	return bwt;
}

int mb_bwt_save(const char *fn, const mb_bwt_t *bwt)
{
	FILE *fp;
	fp = fopen(fn, "wb");
	if (fp == 0) return -1;
	fwrite(MB_MAGIC, 1, 4, fp);
	fwrite(&bwt->sa_bit, 4, 1, fp);
	fwrite(&bwt->primary, 8, 1, fp);
	fwrite(&bwt->L2[1], 8, 4, fp);
	fwrite(bwt->data, 8, bwt->data_len, fp);
	fwrite(&bwt->n_sa, 8, 1, fp);
	if (bwt->sa_bit != (uint32_t)-1 && bwt->n_sa > 0 && bwt->sa)
		fwrite(bwt->sa, 8, bwt->n_sa, fp);
	fclose(fp);
	return 0;
}

mb_bwt_t *mb_bwt_load(const char *fn)
{
	return mb_bwt_load_mt(fn, 1, 0);
}

mb_bwt_t *mb_bwt_load_mt(const char *fn, int n_threads, kom_pfor_f pfor)
{
	FILE *fp;
	char magic[4];
	uint64_t x[5];
	mb_bwt_t *bwt;

	fp = fopen(fn, "rb");
	if (fp == 0) return 0;
	fread(magic, 1, 4, fp);
	if (strncmp(magic, MB_MAGIC, 4) != 0) {
		fclose(fp);
		return 0;
	}
	bwt = mb_bwt_init();
	fread(&bwt->sa_bit, 4, 1, fp);
	fread(x, 8, 5, fp);
	bwt->primary = x[0];
	memcpy(&bwt->L2[1], &x[1], 32);
	bwt->seq_len = bwt->L2[4];
	bwt->data_len = mb_bwt_data_len(bwt->seq_len);
	bwt->data = kom_calloc(uint64_t, bwt->data_len);
	if (kom_fread_par(fp, bwt->data, bwt->data_len << 3, n_threads, pfor) < 0)
		read_huge(fp, bwt->data_len << 3, bwt->data);
	fread(&bwt->n_sa, 8, 1, fp);
	if (bwt->sa_bit != (uint32_t)-1 && bwt->n_sa > 0) {
		uint64_t expected_n_sa = (bwt->seq_len + (1ULL << bwt->sa_bit)) >> bwt->sa_bit;
		if (bwt->n_sa != expected_n_sa) {
			mb_bwt_destroy(bwt);
			return NULL;
		}
		bwt->sa = kom_malloc(uint64_t, bwt->n_sa);
		if (kom_fread_par(fp, bwt->sa, bwt->n_sa << 3, n_threads, pfor) < 0)
			fread(bwt->sa, 8, bwt->n_sa, fp);
	}
	fclose(fp);
	return bwt;
}

mb_bwt_t *mb_bwt_load_mmap(const char *fn, int preload)
{
	uint8_t *base;
	size_t map_len;
	uint64_t data_off, sa_off, n_sa_off, min_len;
	mb_bwt_t *bwt;

	base = (uint8_t*)kom_mmap_file(fn, &map_len, preload);
	if (base == 0) return 0;
	if (map_len < 48 || strncmp((const char*)base, MB_MAGIC, 4) != 0) { kom_munmap(base, map_len); return 0; }

	bwt = mb_bwt_init();
	bwt->mmap = base;
	bwt->mmap_len = map_len;
	bwt->sa_bit = *(const uint32_t*)(base + 4);
	bwt->primary = *(const uint64_t*)(base + 8);
	memcpy(&bwt->L2[1], base + 16, 32);
	bwt->seq_len = bwt->L2[4];
	bwt->data_len = mb_bwt_data_len(bwt->seq_len);

	data_off = 48;
	n_sa_off = data_off + bwt->data_len * 8;
	sa_off = n_sa_off + 8;
	min_len = sa_off; // file must at least hold up to n_sa
	if (map_len < min_len) { mb_bwt_destroy(bwt); return 0; }
	bwt->data = (uint64_t*)(base + data_off);
	bwt->n_sa = *(const uint64_t*)(base + n_sa_off);
	if (bwt->sa_bit != (uint32_t)-1 && bwt->n_sa > 0) {
		uint64_t expected_n_sa = (bwt->seq_len + (1ULL << bwt->sa_bit)) >> bwt->sa_bit;
		if (bwt->n_sa != expected_n_sa || map_len < sa_off + bwt->n_sa * 8) {
			mb_bwt_destroy(bwt);
			return 0;
		}
		bwt->sa = (uint64_t*)(base + sa_off);
	}
	return bwt;
}
