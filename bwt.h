#ifndef MB_BWT_H
#define MB_BWT_H

#include <stdint.h>
#include <stddef.h>

#define MB_MAGIC "MBW\2"

/* A value in mb_bwt_sa_batch()'s input that is already a text position, not an
 * SA row; also set in mb_sai_t::x[1] of a size-1 interval from
 * mb_bwt_smem_batch_pac(), whose x[1] is then the text position of x[0] */
#define MB_SA_KNOWN (1ULL<<63)

typedef struct {
	uint64_t x[2];
	uint64_t size, info;
} mb_sai_t;

typedef struct {
	uint64_t primary; // S^{-1}(0), or the primary index of BWT
	uint64_t L2[5]; // C(), cumulative count
	uint64_t seq_len; // sequence length
	uint64_t data_len;
	uint64_t *data; // BWT
	uint32_t cnt_table[256];
	uint32_t pre_len;
	mb_sai_t *pre;
	// suffix array
	uint32_t sa_bit; // sample rate: 1/(1<<sa_bit)
	uint64_t n_sa;
	uint64_t *sa;
	// mmap
	void *mmap; // mmap pointer
	size_t mmap_len; // mmap length (for munmap)
	// a longer k-mer table beside pre (mb_bwt_cache2): pre2_len-mers, pre_len < pre2_len <= pre_len + 2
	uint32_t pre2_len;
	mb_sai_t *pre2;
} mb_bwt_t;

typedef struct { size_t n, m; mb_sai_t *a; } mb_sai_v;

typedef struct {
	// input before calling
	int32_t min_len, min_occ;
	int32_t st, en;
	const uint8_t *q;
	mb_sai_v *v; // output; v->n is not set to 0
	// internal state
	int32_t stage;
	int32_t x, i;
	uint32_t kmer; // for k-mer cache
	mb_sai_t p;
} mb_smem_entry_t;

/* Fills the slot s with the next job (min_len, min_occ, st, en, q and v) and returns 1, or
 * returns 0 when no job is pending. done is 1 when s holds a job that has just finished, so
 * that its output can generate further jobs, and 0 when s is an idle slot. */
typedef int (*mb_smem_feed_f)(void *data, mb_smem_entry_t *s, int done);

#ifdef __cplusplus
extern "C" {
#endif

extern int mb_verbose;

mb_bwt_t *mb_bwt_init(void);
void mb_bwt_destroy(mb_bwt_t *bwt);
mb_bwt_t *mb_bwt_load_raw(const char *fn); // from raw bwt_gen.c output
int mb_bwt_save(const char *fn, const mb_bwt_t *bwt);
mb_bwt_t *mb_bwt_load(const char *fn);
mb_bwt_t *mb_bwt_load_mt(const char *fn, int n_threads, void (*pfor)(int, void (*)(void*, long, int), void*, long)); // pfor: kt_for() or NULL
mb_bwt_t *mb_bwt_load_mmap(const char *fn, int preload);
mb_bwt_t *mb_bwt_init_from_raw(int is_byte, const void *raw_, uint64_t len, uint64_t primary);

void mb_bwt_cache(mb_bwt_t *bwt, int32_t len);
void mb_bwt_cache2(mb_bwt_t *bwt, int32_t len, int n_threads, void (*pfor)(int, void (*)(void*, long, int), void*, long)); // after mb_bwt_cache(); pfor: kt_for() or NULL
uint64_t mb_bwt_rank11(const mb_bwt_t *bwt, uint64_t k, uint8_t c);
void mb_bwt_rank1a(const mb_bwt_t *bwt, uint64_t k, uint64_t cnt[4]);
void mb_bwt_rank2a(const mb_bwt_t *bwt, uint64_t k, uint64_t l, uint64_t cntk[4], uint64_t cntl[4]);

void mb_bwt_extend(const mb_bwt_t *bwt, const mb_sai_t *ik, mb_sai_t ok[4], int is_back);
void mb_bwt_count_kmer(const mb_bwt_t *bwt, int32_t depth, mb_sai_t *a);
int64_t mb_bwt_smem(const mb_bwt_t *f, uint32_t len, const uint8_t *q, int64_t x, int64_t min_len, int64_t min_occ, mb_sai_t *p);
void mb_bwt_smem_batch(void *km, const mb_bwt_t *bwt, int32_t n, mb_smem_entry_t *a);
void mb_bwt_smem_batch_pac(void *km, const mb_bwt_t *bwt, const uint64_t *pac, int32_t n, mb_smem_entry_t *a);
void mb_bwt_smem_feed(void *km, const mb_bwt_t *bwt, const uint64_t *pac, int32_t n, mb_smem_entry_t *a, mb_smem_feed_f feed, void *data); // feed: NULL or the job source of n slots

void mb_bwt_gen_sa(mb_bwt_t *bwt, uint32_t sa_bit);
uint64_t mb_bwt_sa(const mb_bwt_t *bwt, uint64_t k);
void mb_bwt_sa_batch(void *km, const mb_bwt_t *bwt, int64_t n, uint64_t *x);

static inline void mb_bwt_set_intv(const mb_bwt_t *bwt, int c, mb_sai_t *ik)
{
	ik->x[0] = bwt->L2[c] + 1; // +1 for the missing sentinel
	ik->x[1] = bwt->L2[3-c] + 1;
	ik->size = bwt->L2[c+1] - bwt->L2[c];
	ik->info = 0;
}
#ifdef __cplusplus
}
#endif

#endif
