#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdarg.h>
#include <pthread.h>
#include <sys/stat.h>
#define __STDC_LIMIT_MACROS
#include "bseq.h"
#include "kseq.h"

static int bseq_gzread(mb_bseq_file_t *fp, void *buf, int len);
KSEQ_INIT2(, mb_bseq_file_t*, bseq_gzread)

#define kvec_t(type) struct { size_t n, m; type *a; }

#define kv_resize(type, v, s) do { \
		if ((v).m < (s)) { \
			(v).m = (s); \
			(v).a = (type*)realloc((v).a, sizeof(type) * (v).m); \
		} \
	} while (0)

#define kv_push(type, v, x) do { \
		if ((v).n == (v).m) { \
			(v).m += ((v).m>>1) + 16; \
			(v).a = (type*)realloc((v).a, sizeof(type) * (v).m); \
		} \
		(v).a[(v).n++] = (x); \
	} while (0)

#define kv_pushp(type, v, p) do { \
		if ((v).n == (v).m) { \
			(v).m += ((v).m>>1) + 16; \
			(v).a = (type*)realloc((v).a, sizeof(type) * (v).m); \
		} \
		*(p) = &(v).a[(v).n++]; \
	} while (0)

#define CHECK_PAIR_THRES 1000000

static kstring_t *bseq_warn_buf = 0; // when set, warnings are appended here instead of printed

void mb_bseq_warn_buffer(kstring_t *buf)
{
	bseq_warn_buf = buf;
}

static void (*bseq_hook)(void *data, const mb_bseq1_t *a, int n) = 0; // see mb_bseq_read_hook()
static void *bseq_hook_data = 0;

void mb_bseq_read_hook(void (*hook)(void *data, const mb_bseq1_t *a, int n), void *data)
{
	bseq_hook = hook, bseq_hook_data = data;
}

static void bseq_warn(const char *fmt, ...)
{
	va_list ap;
	va_start(ap, fmt);
	if (bseq_warn_buf) {
		kstring_t *s = bseq_warn_buf;
		va_list aq;
		int len;
		va_copy(aq, ap);
		len = vsnprintf(0, 0, fmt, aq);
		va_end(aq);
		if (len > 0) {
			if (s->l + len + 1 > s->m) {
				s->m = s->l + len + 1;
				s->m += s->m >> 1;
				s->s = (char*)realloc(s->s, s->m);
			}
			vsnprintf(s->s + s->l, len + 1, fmt, ap);
			s->l += len;
		}
	} else vfprintf(stderr, fmt, ap);
	va_end(ap);
}

/*
 * Inflating ahead of the parser
 *
 * mb_bseq_ahead(fp, 1) starts a thread that makes the gzread() calls kseq
 * would make on fp, one per 16384-byte chunk (kseq's buffer size), in the same
 * order, and queues each chunk with gzread()'s return value; kseq then takes
 * the chunks from the queue in order (bseq_gzread()), so it sees the same
 * bytes and return values. mb_bseq_ahead(fp, 0) stops the thread after the
 * block it is filling; the chunks it has queued are read first, then
 * bseq_gzread() calls gzread() itself again. The thread is started only for a
 * regular file, which gzread() never waits on.
 */
#define BSEQ_CHUNK      16384 // kseq's buffer size in KSEQ_INIT2()
#define BSEQ_BLK_CHUNKS 16    // chunks handed over at a time
#define BSEQ_N_BLK      8     // blocks queued at most

typedef struct {
	int32_t n, ret[BSEQ_BLK_CHUNKS]; // number of chunks and gzread()'s return value for each
	unsigned char *data;
} bseq_blk_t;

typedef struct {
	pthread_t tid;
	pthread_mutex_t lock;
	pthread_cond_t cv;
	int32_t stop, done;   // done: the thread has queued its last block
	int32_t head, cnt;    // queued blocks: blk[head], ..., blk[(head+cnt-1)%BSEQ_N_BLK]
	int32_t pos;          // the next chunk of blk[head]
	bseq_blk_t blk[BSEQ_N_BLK];
	unsigned char *buf;   // gzread()'s buffer, reused for every call as kseq reuses its own
} bseq_ahead_t;

struct mb_bseq_file_s {
	gzFile fp;
	kseq_t *ks;
	mb_bseq1_t s;
	int32_t is_reg;
	bseq_ahead_t *ah;
};

static void *bseq_ahead_worker(void *data)
{
	mb_bseq_file_t *fp = (mb_bseq_file_t*)data;
	bseq_ahead_t *ah = fp->ah;
	bseq_blk_t *b = 0;
	int32_t last = 0;
	for (;;) {
		pthread_mutex_lock(&ah->lock);
		if (b) ++ah->cnt, pthread_cond_signal(&ah->cv);
		while (!ah->stop && !last && ah->cnt == BSEQ_N_BLK)
			pthread_cond_wait(&ah->cv, &ah->lock);
		if (ah->stop || last) {
			ah->done = 1;
			pthread_cond_signal(&ah->cv);
			pthread_mutex_unlock(&ah->lock);
			break;
		}
		b = &ah->blk[(ah->head + ah->cnt) % BSEQ_N_BLK];
		pthread_mutex_unlock(&ah->lock);
		for (b->n = 0; b->n < BSEQ_BLK_CHUNKS && !last; ++b->n) {
			int r;
			r = gzread(fp->fp, ah->buf, BSEQ_CHUNK);
			memcpy(b->data + (size_t)b->n * BSEQ_CHUNK, ah->buf, r < 0? BSEQ_CHUNK : r); // on an error, kseq's buffer holds what gzread() left there
			b->ret[b->n] = r;
			if (r < BSEQ_CHUNK) last = 1; // kseq makes no gzread() call after a short one
		}
	}
	return 0;
}

static void bseq_ahead_free(bseq_ahead_t *ah)
{
	int32_t i;
	pthread_mutex_destroy(&ah->lock);
	pthread_cond_destroy(&ah->cv);
	for (i = 0; i < BSEQ_N_BLK; ++i) free(ah->blk[i].data);
	free(ah->buf);
	free(ah);
}

static void bseq_ahead_stop(mb_bseq_file_t *fp)
{
	bseq_ahead_t *ah = fp->ah;
	if (ah == 0 || ah->stop) return;
	pthread_mutex_lock(&ah->lock);
	ah->stop = 1;
	pthread_cond_signal(&ah->cv);
	pthread_mutex_unlock(&ah->lock);
	pthread_join(ah->tid, 0);
	if (ah->cnt == 0) bseq_ahead_free(ah), fp->ah = 0;
}

void mb_bseq_ahead(mb_bseq_file_t *fp, int on)
{
	bseq_ahead_t *ah;
	int32_t i;
	if (!on) {
		bseq_ahead_stop(fp);
		return;
	}
	if (fp->ah || !fp->is_reg) return;
	ah = (bseq_ahead_t*)calloc(1, sizeof(bseq_ahead_t));
	ah->buf = (unsigned char*)malloc(BSEQ_CHUNK);
	for (i = 0; i < BSEQ_N_BLK; ++i)
		ah->blk[i].data = (unsigned char*)malloc((size_t)BSEQ_BLK_CHUNKS * BSEQ_CHUNK);
	pthread_mutex_init(&ah->lock, 0);
	pthread_cond_init(&ah->cv, 0);
	fp->ah = ah;
	if (pthread_create(&ah->tid, 0, bseq_ahead_worker, fp) != 0) {
		bseq_ahead_free(ah);
		fp->ah = 0;
	}
}

static int bseq_gzread(mb_bseq_file_t *fp, void *buf, int len)
{
	bseq_ahead_t *ah = fp->ah;
	bseq_blk_t *b;
	int r;
	if (ah == 0) return gzread(fp->fp, buf, len);
	assert(len == BSEQ_CHUNK);
	pthread_mutex_lock(&ah->lock);
	while (ah->cnt == 0 && !ah->done)
		pthread_cond_wait(&ah->cv, &ah->lock);
	if (ah->cnt == 0) { // the thread was stopped and everything it read has been taken
		pthread_mutex_unlock(&ah->lock);
		if (!ah->stop) pthread_join(ah->tid, 0); // it stopped at a short read; kseq would not call this again
		bseq_ahead_free(ah);
		fp->ah = 0;
		return gzread(fp->fp, buf, len);
	}
	b = &ah->blk[ah->head];
	pthread_mutex_unlock(&ah->lock);
	r = b->ret[ah->pos];
	memcpy(buf, b->data + (size_t)ah->pos * BSEQ_CHUNK, r < 0? BSEQ_CHUNK : r);
	if (++ah->pos == b->n) { // hand the block back
		ah->pos = 0;
		pthread_mutex_lock(&ah->lock);
		ah->head = (ah->head + 1) % BSEQ_N_BLK;
		--ah->cnt;
		pthread_cond_signal(&ah->cv);
		pthread_mutex_unlock(&ah->lock);
	}
	return r;
}

mb_bseq_file_t *mb_bseq_open(const char *fn)
{
	mb_bseq_file_t *fp;
	gzFile f;
	struct stat st;
	int is_stdin = !(fn && strcmp(fn, "-"));
	f = !is_stdin? gzopen(fn, "r") : gzdopen(0, "r");
	if (f == 0) return 0;
	fp = (mb_bseq_file_t*)calloc(1, sizeof(mb_bseq_file_t));
	fp->fp = f;
	fp->is_reg = (is_stdin? fstat(0, &st) : stat(fn, &st)) == 0 && S_ISREG(st.st_mode);
	fp->ks = kseq_init(fp);
	return fp;
}

void mb_bseq_close(mb_bseq_file_t *fp)
{
	if (fp->ah) {
		bseq_ahead_stop(fp);
		if (fp->ah) bseq_ahead_free(fp->ah);
	}
	kseq_destroy(fp->ks);
	gzclose(fp->fp);
	free(fp);
}

static inline char *kstrdup(const kstring_t *s)
{
	char *t;
	t = (char*)malloc(s->l + 1);
	memcpy(t, s->s, s->l + 1);
	return t;
}

static inline void kseq2bseq(kseq_t *ks, mb_bseq1_t *s, int with_qual, int with_comment)
{
	int i;
	if (ks->name.l == 0)
		bseq_warn("[WARNING]\033[1;31m empty sequence name in the input.\033[0m\n");
	s->name = kstrdup(&ks->name);
	s->seq = kstrdup(&ks->seq);
	for (i = 0; i < (int)ks->seq.l; ++i) // convert U to T
		if (s->seq[i] == 'u' || s->seq[i] == 'U')
			--s->seq[i];
	s->qual = with_qual && ks->qual.l? kstrdup(&ks->qual) : 0;
	s->comment = with_comment && ks->comment.l? kstrdup(&ks->comment) : 0;
	s->l_seq = ks->seq.l;
}

mb_bseq1_t *mb_bseq_read(mb_bseq_file_t *fp, int64_t chunk_size, int with_qual, int with_comment, int frag_mode, int min_cnt, int64_t max_chunk_size, int *n_)
{
	int64_t size = 0;
	int ret;
	kvec_t(mb_bseq1_t) a = {0,0,0};
	kseq_t *ks = fp->ks;
	*n_ = 0;
	if (fp->s.seq) {
		kv_resize(mb_bseq1_t, a, 256);
		kv_push(mb_bseq1_t, a, fp->s);
		size = fp->s.l_seq;
		memset(&fp->s, 0, sizeof(mb_bseq1_t));
		if (bseq_hook) bseq_hook(bseq_hook_data, a.a, a.n);
	}
	if (max_chunk_size < chunk_size)
		max_chunk_size = chunk_size;
	while ((ret = kseq_read(ks)) >= 0) {
		int32_t to_stop = 0;
		mb_bseq1_t *s;
		assert(ks->seq.l <= INT32_MAX);
		if (a.m == 0) kv_resize(mb_bseq1_t, a, 256);
		kv_pushp(mb_bseq1_t, a, &s);
		kseq2bseq(ks, s, with_qual, with_comment);
		if (bseq_hook) bseq_hook(bseq_hook_data, a.a, a.n);
		size += ks->seq.l;
		if (chunk_size <= 0 || max_chunk_size <= 0) to_stop = 1;
		else if (size >= max_chunk_size) to_stop = 1;
		else if (size >= chunk_size && a.n >= min_cnt) to_stop = 1;
		if (to_stop) {
			if (frag_mode && a.a[a.n-1].l_seq < CHECK_PAIR_THRES) {
				while ((ret = kseq_read(ks)) >= 0) {
					kseq2bseq(ks, &fp->s, with_qual, with_comment);
					size += ks->seq.l;
					if (mb_qname_same(fp->s.name, a.a[a.n-1].name)) {
						kv_push(mb_bseq1_t, a, fp->s);
						memset(&fp->s, 0, sizeof(mb_bseq1_t));
						if (bseq_hook) bseq_hook(bseq_hook_data, a.a, a.n);
					} else break;
				}
			}
			break;
		}
	}
	if (ret < -1) {
		if (a.n) bseq_warn("[WARNING]\033[1;31m failed to parse the FASTA/FASTQ record next to '%s'. Continue anyway.\033[0m\n", a.a[a.n-1].name);
		else bseq_warn("[WARNING]\033[1;31m failed to parse the first FASTA/FASTQ record. Continue anyway.\033[0m\n");
	}
	*n_ = a.n;
	return a.a;
}

mb_bseq1_t *mb_bseq_read_frag(int n_fp, mb_bseq_file_t **fp, int64_t chunk_size, int with_qual, int with_comment, int *n_)
{
	int i;
	int64_t size = 0;
	kvec_t(mb_bseq1_t) a = {0,0,0};
	*n_ = 0;
	if (n_fp < 1) return 0;
	while (1) {
		int n_read = 0;
		for (i = 0; i < n_fp; ++i)
			if (kseq_read(fp[i]->ks) >= 0)
				++n_read;
		if (n_read < n_fp) {
			if (n_read > 0)
				bseq_warn("[W::%s]\033[1;31m query files have different number of records; extra records skipped.\033[0m\n", __func__);
			break; // some file reaches the end
		}
		if (a.m == 0) kv_resize(mb_bseq1_t, a, 256);
		for (i = 0; i < n_fp; ++i) {
			mb_bseq1_t *s;
			kv_pushp(mb_bseq1_t, a, &s);
			kseq2bseq(fp[i]->ks, s, with_qual, with_comment);
			size += s->l_seq;
		}
		if (bseq_hook) bseq_hook(bseq_hook_data, a.a, a.n);
		if (size >= chunk_size) break;
	}
	*n_ = a.n;
	return a.a;
}

int mb_bseq_eof(mb_bseq_file_t *fp)
{
	return (ks_eof(fp->ks->f) && fp->s.seq == 0);
}
