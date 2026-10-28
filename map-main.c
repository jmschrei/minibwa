#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <errno.h>
#include <zlib.h>
#include "kommon.h"
#include "mbpriv.h"
#include "bseq.h"
#include "kalloc.h"
#include "kthread.h"
#include "ketopt.h"
#include "kseq.h"
KSTREAM_INIT(gzFile, gzread, 0x10000)

typedef struct {
	int32_t n_fp, n_threads;
	int64_t n_base, n_seq, n_batch, n_batch_eof; // n_batch_eof: batches read when step 0 met the end of input; 0 before
	int64_t mb_size; // mini-batch size
	const mb_opt_t *opt;
	mb_bseq_file_t **fp;
	const mb_idx_t *idx;
	FILE *fp_out;
	mb_bseq1_t *pre_seq; // the first batch, read while the index loaded (mb_prefetch_t)
	int32_t pre, pre_n; // pre: step 0 has yet to take pre_seq[0..pre_n-1] as its first batch
	void *early; // sub-batches of the first batch mapped while it was read (mb_early_t), until step 1 takes them; 0 if none
} pipeline_t;

typedef struct {
	const pipeline_t *p;
    int32_t n_seq, n_frag, n_sb, n_pe, n_out, out_done; // out_done: out[] formatted in step 1
	int32_t *n_hit, *seg_off, *seg_cnt;
	int32_t *sb_off, *sb_cnt;
	mb_pestat_t pes[4];
	mb_bseq1_t *seq;
	mb_hit_t **hit;
	mb_tbuf_t **tbuf;
	kstring_t *out; // out[n_out]: SAM/PAF text of chunks of fragments
	int64_t batch; // the batch's number, from 0
} step_t;

static void worker_for_se_batch(void *data, long i, int tid)
{
	step_t *s = (step_t*)data;
	const mb_opt_t *opt = s->p->opt;
	const mb_idx_t *idx = s->p->idx;
	mb_tbuf_t *b = s->tbuf[tid];
	void *km;
	int64_t tot;
	int32_t n, j, k, l, p, *len;
	uint8_t **seq, *buf;
	mb_sai_v *sai;

	km = mb_tbuf_km(b);
	for (k = 0, n = 0, tot = 0; k < s->sb_cnt[i]; ++k) {
		n += s->seg_cnt[s->sb_off[i] + k];
		int32_t off = s->seg_off[s->sb_off[i] + k];
		int32_t cnt = s->seg_cnt[s->sb_off[i] + k];
		for (j = 0; j < cnt; ++j)
			tot += s->seq[off + j].l_seq;
	}
	buf = Kmalloc(km, uint8_t, tot);
	len = Kmalloc(km, int32_t, n);
	seq = Kmalloc(km, uint8_t*, n);
	sai = Kcalloc(km, mb_sai_v, n);
	for (k = p = 0, tot = 0; k < s->sb_cnt[i]; ++k) {
		int32_t off = s->seg_off[s->sb_off[i] + k];
		int32_t cnt = s->seg_cnt[s->sb_off[i] + k];
		for (j = 0; j < cnt; ++j) {
			const mb_bseq1_t *t = &s->seq[off + j];
			l2b_meth_t mt = !idx->is_meth? L2B_METH_NONE : (j&1) == 0? L2B_METH_C2T : L2B_METH_G2A;
			len[p] = t->l_seq;
			seq[p] = &buf[tot], tot += t->l_seq;
			for (l = 0; l < t->l_seq; ++l)
				seq[p][l] = kom_nt4_table[(uint8_t)t->seq[l]];
			if (mt != L2B_METH_NONE) l2b_meth_convert(mt, len[p], seq[p]);
			++p;
		}
	}
	assert(p == n);
	uint64_t *sa = mb_seed_sa_batch(km, idx, n, len, seq, opt->min_len, opt->max_sub_occ, opt->max_occ, sai), *sa_p = sa; // mb_opt_adap() keeps max_occ
	kfree(km, seq);
	kfree(km, len);
	kfree(km, buf);

	for (k = p = 0; k < s->sb_cnt[i]; ++k) {
		int32_t off = s->seg_off[s->sb_off[i] + k];
		int32_t cnt = s->seg_cnt[s->sb_off[i] + k];
		for (j = 0; j < cnt; ++j, ++p) {
			const mb_bseq1_t *t = &s->seq[off + j];
			mb_opt_t opt_adap;
			l2b_meth_t mt = !idx->is_meth? L2B_METH_NONE : (j&1) == 0? L2B_METH_C2T : L2B_METH_G2A;
			mb_opt_adap(opt, t->l_seq, &opt_adap);
			s->hit[off+j] = mb_map_sai_sa(&opt_adap, idx, t->l_seq, t->seq, mt, &sai[p], sa? &sa_p : 0, &s->n_hit[off+j], b, t->name);
		}
	}
	kfree(km, sa);
	kfree(km, sai);
	mb_tbuf_reset(b, opt->cap_kalloc);
}

static void worker_for_pe(void *data, long i, int tid)
{
	step_t *s = (step_t*)data;
	mb_tbuf_t *b = s->tbuf[tid];
	int32_t off, r, len[2];
	char *seq[2];
	if (s->seg_cnt[i] != 2) return;
	off = s->seg_off[i];
	if (kom_dbg_flag & MB_DBG_QNAME) fprintf(stderr, "QP\t%s\t%d\n", s->seq[off].name, tid);
	for (r = 0; r < 2; ++r) {
		seq[r] = s->seq[off + r].seq;
		len[r] = s->seq[off + r].l_seq;
	}
	mb_pair(mb_tbuf_km(b), s->p->opt, s->p->idx->l2b, &s->n_hit[off], &s->hit[off], s->pes, len, seq);
}

static void worker_for_out(void *data, long c, int tid); // formats chunk c of a batch's output; below, after worker_pipeline()
static int32_t mb_early_merge(pipeline_t *p, step_t *s); // takes the hits of the first batch mapped early; below, with mb_early_t

static void *worker_pipeline(void *shared, int step, void *in)
{
	const int min_read_cnt = 40000; // should be smaller than opt->mb_size / (read_len * 2) for typical paired-end reads
	int i, j, k;
    pipeline_t *p = (pipeline_t*)shared;
	const mb_opt_t *opt = p->opt;
    if (step == 0) { // step 0: read sequences
		int with_qual = !(opt->flag & MB_F_PAF);
		int with_comment = !!(opt->flag & MB_F_COPY_COMMENT);
		int frag_mode = (p->n_fp > 1 || !!(opt->flag & MB_F_PE));
        step_t *s;
        s = kom_calloc(step_t, 1);
		if (p->pre) { // the first batch, read while the index loaded
			s->seq = p->pre_seq;
			s->n_seq = p->pre_n;
			p->pre = 0;
		} else if (p->n_fp > 1) s->seq = mb_bseq_read_frag(p->n_fp, p->fp, p->mb_size, with_qual, with_comment, &s->n_seq);
		else s->seq = mb_bseq_read(p->fp[0], p->mb_size, with_qual, with_comment, frag_mode, min_read_cnt, opt->max_mb_size, &s->n_seq);
		if (s->seq) {
			int32_t sb_len, sb_off;
			s->p = p, s->batch = p->n_batch++;
			for (i = 0; i < 4; ++i) s->pes[i].failed = 1;
			for (i = 0; i < s->n_seq; ++i)
				s->seq[i].id = p->n_seq++;
			s->tbuf = kom_calloc(mb_tbuf_t*, opt->n_thread);
			for (i = 0; i < opt->n_thread; ++i)
				s->tbuf[i] = mb_tbuf_init(opt->flag&MB_F_NO_KALLOC);
			s->n_hit = kom_calloc(int32_t, 5 * s->n_seq); // maybe over allocation as the following 4 members may not need this much
			s->seg_off = s->n_hit   + s->n_seq;
			s->seg_cnt = s->seg_off + s->n_seq;
			s->sb_off  = s->seg_cnt + s->n_seq;
			s->sb_cnt  = s->sb_off  + s->n_seq;
			s->hit = kom_calloc(mb_hit_t*, s->n_seq);
			// set seg_cnt[] and seg_off[]
			for (i = 1, j = 0; i <= s->n_seq; ++i) {
				if (i == s->n_seq || !frag_mode || !mb_qname_same(s->seq[i-1].name, s->seq[i].name)) {
					assert(i - j <= 2);
					s->seg_cnt[s->n_frag] = i - j;
					s->seg_off[s->n_frag++] = j;
					if (i - j == 2) s->n_pe++;
					j = i;
				}
			}
			if (s->n_pe > 0) { // trim trailing "/[12]" for paired-end reads
				for (i = 0; i < s->n_frag; ++i) {
					int32_t j0, j1, l0, l1;
					if (s->seg_cnt[i] != 2) continue;
					j0 = s->seg_off[i] + 0, l0 = strlen(s->seq[j0].name);
					j1 = s->seg_off[i] + 1, l1 = strlen(s->seq[j1].name);
					if (l0 < 3 || l0 != l1) continue;
					if (s->seq[j0].name[l0 - 1] != s->seq[j1].name[l1 - 1] && s->seq[j0].name[l0 - 2] == '/')
						s->seq[j0].name[l0 - 2] = s->seq[j1].name[l1 - 2] = 0; // truncate
				}
			}
			// set sb_cnt[] and sb_off[]
			for (i = 0, sb_len = sb_off = 0; i < s->n_frag; ++i) {
				if (sb_len >= opt->sb_len || i - sb_off >= opt->sb_seq) {
					s->sb_off[s->n_sb] = sb_off;
					s->sb_cnt[s->n_sb++] = i - sb_off;
					sb_len = 0;
					sb_off = i;
				}
				for (j = 0; j < s->seg_cnt[i]; ++j)
					sb_len += s->seq[s->seg_off[i] + j].l_seq;
			}
			s->sb_off[s->n_sb] = sb_off;
			s->sb_cnt[s->n_sb++] = i - sb_off;
			return s;
		} else {
			__sync_bool_compare_and_swap(&p->n_batch_eof, 0, p->n_batch); // no batch follows batch n_batch-1
			free(s);
		}
    } else if (step == 1) { // step 1: map
		step_t *s = (step_t*)in;
		int32_t n_sb = s->n_sb;
		if (p->early) // the first batch, of which some sub-batches were mapped while it was read
			n_sb = mb_early_merge(p, s);
		kt_for(opt->n_thread, worker_for_se_batch, in, n_sb);
		if ((opt->flag & MB_F_PE) && s->n_frag < s->n_seq && !(opt->flag & (MB_F_NO_PAIRING|MB_F_NO_ALN))) { // PE mode
			if ((opt->flag & MB_F_PE_PREDEF) || s->n_pe < 20) { // use predefined PE stats
				s->pes[1].failed = 0;
				s->pes[1].avg = opt->pe_avg, s->pes[1].std = opt->pe_std;
				s->pes[1].lo = opt->pe_lo, s->pes[1].hi = opt->pe_hi;
			} else { // estimate PE stats from data
				void *km;
				km = opt->flag & MB_F_NO_KALLOC? 0 : km_init();
				mb_pestat(km, opt, s->n_frag, s->seg_off, s->seg_cnt, s->n_hit, s->hit, s->pes);
				if (km) km_destroy(km);
			}
			kt_for(opt->n_thread, worker_for_pe, in, s->n_frag);
		}
		s->n_out = (s->n_frag + 255) / 256; // chunks of about 256 fragments
		s->out = kom_calloc(kstring_t, s->n_out);
		if (s->batch + 1 == __sync_fetch_and_add(&p->n_batch_eof, 0)) { // the last batch: no mapping follows that step 2 could overlap, so format on all threads here
			kt_for(opt->n_thread, worker_for_out, in, s->n_out);
			s->out_done = 1;
		}
		return in;
    } else if (step == 2) { // step 2: output
		step_t *s = (step_t*)in;
		kstring_t buf = {0,0,0};
		int64_t tot_len = 0;
		for (k = 0; k < s->n_out; ++k) {
			if (!s->out_done) { // not formatted in step 1: format chunk k here, reusing one buffer
				buf.l = 0, s->out[k] = buf;
				worker_for_out(s, k, 0);
				buf = s->out[k];
			}
			if (s->out[k].l > 0) fwrite(s->out[k].s, 1, s->out[k].l, s->p->fp_out);
			if (s->out_done) free(s->out[k].s);
		}
		free(buf.s);
		for (i = 0; i < s->n_seq; ++i)
			tot_len += s->seq[i].l_seq;
		for (i = 0; i < opt->n_thread; ++i)
			mb_tbuf_destroy(s->tbuf[i]);
		free(s->tbuf);
		free(s->out); free(s->hit); free(s->n_hit); free(s->seq);
		if (kom_verbose >= 3)
			fprintf(stderr, "[M::%s::%.3f*%.2f] mapped %ld bp in %ld sequences\n", __func__, kom_realtime(), kom_percent_cpu(), (long)tot_len, (long)s->n_seq);
		free(s);
	}
    return 0;
}

static void worker_for_out(void *data, long c, int tid) // formats chunk c of the fragments into s->out[c], then frees their reads and hits
{
	step_t *s = (step_t*)data;
	const mb_opt_t *opt = s->p->opt;
	const l2b_t *l2b = s->p->idx->l2b;
	void *km = mb_tbuf_km(s->tbuf[tid]);
	kstring_t *out = &s->out[c];
	int32_t i, j, k, k_st = (int64_t)c * s->n_frag / s->n_out, k_en = (int64_t)(c + 1) * s->n_frag / s->n_out;

	for (k = k_st; k < k_en; ++k) {
		int32_t seg_st = s->seg_off[k], seg_en = s->seg_off[k] + s->seg_cnt[k];
		for (i = seg_st; i < seg_en; ++i) {
			mb_bseq1_t *t = &s->seq[i];
			int32_t mate_qlen = 0; // mate's l_seq for MC:Z; 0 suppresses MC/MQ
			if (seg_en - seg_st > 1) {
				int32_t mate_idx = i != seg_en - 1? i + 1 : seg_st; // wrap around if i is the last segment
				mate_qlen = s->seq[mate_idx].l_seq;
			}
			if (s->n_hit[i] > 0) { // the query has at least one hit
				int32_t n_sec = 0;
				for (j = 0; j < s->n_hit[i]; ++j) {
					const mb_hit_t *h = &s->hit[i][j];
					if (h->parent == h->id || n_sec < opt->out_n) {
						if (h->parent != h->id) {
							const mb_hit_t *p = &s->hit[i][h->parent];
							if (p->p && h->p) {
								if (h->p->dp_max < (double)opt->out_s * p->p->dp_max) continue;
							} else {
								if (h->score < (double)opt->out_s * p->score) continue;
							}
						}
						mb_format(km, out, l2b, t, seg_en - seg_st, &s->n_hit[seg_st], &s->hit[seg_st], j, opt, i - seg_st, mate_qlen);
						n_sec += (h->parent != h->id);
					}
				}
			} else if (!(opt->flag & MB_F_NO_UNMAP)) {
				mb_format(km, out, l2b, t, seg_en - seg_st, &s->n_hit[seg_st], &s->hit[seg_st], -1, opt, i - seg_st, mate_qlen);
			}
		}
		for (i = seg_st; i < seg_en; ++i) {
			for (j = 0; j < s->n_hit[i]; ++j) free(s->hit[i][j].p);
			free(s->hit[i]);
			free(s->seq[i].seq); free(s->seq[i].name);
			if (s->seq[i].qual) free(s->seq[i].qual);
			if (s->seq[i].comment) free(s->seq[i].comment);
		}
	}
}

static mb_bseq_file_t **mb_open_bseqs(int n, const char **fn)
{
	mb_bseq_file_t **fp;
	int32_t i, j;
	fp = kom_calloc(mb_bseq_file_t*, n);
	for (i = 0; i < n; ++i) {
		if ((fp[i] = mb_bseq_open(fn[i])) == 0) {
			if (kom_verbose >= 1)
				fprintf(stderr, "ERROR: failed to open file '%s': %s\n", fn[i], strerror(errno));
			for (j = 0; j < i; ++j)
				mb_bseq_close(fp[j]);
			free(fp);
			return 0;
		}
	}
	return fp;
}

/*
 * Reading the first batch while the index loads
 *
 * main_map() and main_mem() start a thread that opens the read files and reads
 * the first batch, with the same mb_bseq_read()/mb_bseq_read_frag() call and
 * arguments step 0 of worker_pipeline() uses, before they load the index;
 * mb_map_file_pre() hands that batch to kt_pipeline() as step 0's first
 * result, so batch boundaries do not change. Nothing of it is printed early:
 * a file that fails to open is reported after the SAM header is written, as
 * mb_open_bseqs() reports it, and bseq's warnings are kept and printed where
 * step 0 prints them. On an exit before that (a missing index, a bad -R, an
 * output file that cannot be opened) the thread is not waited for and prints
 * nothing. It is skipped at -t 1, and when an input is missing or is the
 * output file, so that opening the output before the reads still decides.
 * While it reads the first batch, each input that is a regular file is
 * inflated on a thread of its own (mb_bseq_ahead() in bseq.c), so R1 and R2
 * inflate in parallel with each other and with the parsing.
 */
#include <pthread.h>
#include <sys/stat.h>

/*
 * Mapping the first batch while it is read
 *
 * When the index has loaded before the thread below has read the first batch,
 * mb_map_file_pre() starts n_thread-1 workers that map the batch's sub-batches
 * with worker_for_se_batch() as the reader completes them. From then on the
 * reader copies each record it takes, in mb_early_hook(), which bseq calls,
 * into blocks that never move (the batch's own array moves as it grows), and
 * the main thread groups the copies into fragments and the fragments into
 * sub-batches by step 0's rules and queues each sub-batch once the record after
 * it is known. A pair's names are mapped trimmed of "/1" and "/2" as step 0
 * trims them, in a copy, since step 0 groups the untrimmed names again. The
 * workers stop when the read ends; step 1 takes the hits of the sub-batches
 * they mapped after checking each against step 0's fragments and sub-batches,
 * and maps the others; mb_pestat() and pairing then run on the whole batch.
 * It is off with --dbg-* options, whose stderr lines it would move ahead of
 * bseq's warnings.
 */
#define MB_EARLY_BLK_SHIFT 16
#define MB_EARLY_MAX_BLK   (1<<14)

typedef struct {
	int32_t r0, n_rec, n_frag; // records r0..r0+n_rec-1 of the batch, in n_frag fragments
	int32_t done, taken; // done: mapped; taken: its hits are step 1's
	int32_t *seg_off, *seg_cnt, *n_hit; // seg_off[]: from r0
	mb_hit_t **hit;
	mb_bseq1_t *seq; // copies of the records, sharing their names and sequences
} mb_early_job_t;

struct mb_early_s;

typedef struct {
	struct mb_early_s *e;
	int32_t tid;
} mb_early_worker_t;

typedef struct mb_early_s {
	const mb_opt_t *opt;
	const pipeline_t *pl;
	int32_t frag_mode;
	int32_t want, n_pub; // want: the index has loaded; n_pub: records in blk[]; both read and written with __atomic builtins
	// the reader's: records copied, the count at its last wake-up of the main thread, bases copied since then, blk[] is full
	int32_t n_copy, n_sig, full;
	int64_t l_sig;
	mb_bseq1_t **blk; // the copies, 1<<MB_EARLY_BLK_SHIFT records per block
	// the main thread's: records grouped, the open fragment's first record, fragments closed,
	// the open sub-batch's first fragment and length (as step 0's sb_off and sb_len), its fragments
	int32_t n_proc, j0, n_frag, sb_off, sb_len, broken;
	int32_t n_f, m_f, *f_off, *f_cnt;
	// under mtx; cv wakes the workers, cv_pub the main thread
	pthread_mutex_t mtx;
	pthread_cond_t cv, cv_pub;
	int32_t reading, n_job, m_job, next;
	mb_early_job_t **job;
	// the workers
	int32_t n_worker, n_tbuf;
	pthread_t *tid;
	mb_early_worker_t *w;
	mb_tbuf_t **tbuf;
} mb_early_t;

static mb_early_t *mb_early_init(const mb_opt_t *opt, int32_t frag_mode)
{
	mb_early_t *e;
	e = kom_calloc(mb_early_t, 1);
	e->opt = opt, e->frag_mode = frag_mode, e->reading = 1;
	pthread_mutex_init(&e->mtx, 0);
	pthread_cond_init(&e->cv, 0);
	pthread_cond_init(&e->cv_pub, 0);
	return e;
}

static void mb_early_destroy(mb_early_t *e)
{
	int32_t k, r, j;
	for (k = 0; k < e->n_job; ++k) {
		mb_early_job_t *job = e->job[k];
		if (job->done && !job->taken) {
			for (r = 0; r < job->n_rec; ++r) {
				for (j = 0; j < job->n_hit[r]; ++j) free(job->hit[r][j].p);
				free(job->hit[r]);
			}
		}
		free(job->seq); free(job->seg_off); free(job->hit); free(job);
	}
	if (e->blk) {
		for (k = 0; k < MB_EARLY_MAX_BLK && e->blk[k]; ++k)
			free(e->blk[k]);
		free(e->blk);
	}
	free(e->job); free(e->f_off); free(e->f_cnt);
	free(e->tid); free(e->w); free(e->tbuf);
	pthread_mutex_destroy(&e->mtx);
	pthread_cond_destroy(&e->cv);
	pthread_cond_destroy(&e->cv_pub);
	free(e);
}

static void mb_early_hook(void *data, const mb_bseq1_t *a, int n) // the reader: bseq calls it with the batch so far, a[0..n-1]
{
	mb_early_t *e = (mb_early_t*)data;
	const int32_t bs = 1<<MB_EARLY_BLK_SHIFT;
	if (e->full || !__atomic_load_n(&e->want, __ATOMIC_ACQUIRE)) return;
	if (e->blk == 0) e->blk = kom_calloc(mb_bseq1_t*, MB_EARLY_MAX_BLK);
	while (e->n_copy < n) {
		int32_t b = e->n_copy >> MB_EARLY_BLK_SHIFT, o = e->n_copy & (bs - 1), l = n - e->n_copy;
		if (b == MB_EARLY_MAX_BLK) {
			e->full = 1;
			break;
		}
		if (e->blk[b] == 0) e->blk[b] = kom_malloc(mb_bseq1_t, bs);
		if (l > bs - o) l = bs - o;
		memcpy(&e->blk[b][o], &a[e->n_copy], l * sizeof(mb_bseq1_t));
		for (o = 0; o < l; ++o)
			e->l_sig += a[e->n_copy + o].l_seq;
		e->n_copy += l;
	}
	__atomic_store_n(&e->n_pub, e->n_copy, __ATOMIC_RELEASE);
	if (e->n_copy - e->n_sig >= 256 || e->l_sig >= e->opt->sb_len || e->full) { // wake the main thread now and then
		pthread_mutex_lock(&e->mtx);
		pthread_cond_signal(&e->cv_pub);
		pthread_mutex_unlock(&e->mtx);
		e->n_sig = e->n_copy, e->l_sig = 0;
	}
}

static inline const mb_bseq1_t *mb_early_rec(const mb_early_t *e, int32_t i)
{
	return &e->blk[i >> MB_EARLY_BLK_SHIFT][i & ((1<<MB_EARLY_BLK_SHIFT) - 1)];
}

static void mb_early_push(mb_early_t *e) // queues the open sub-batch
{
	mb_early_job_t *job;
	int32_t f, r, r0 = e->f_off[0], n_rec = e->f_off[e->n_f - 1] + e->f_cnt[e->n_f - 1] - r0;
	job = kom_calloc(mb_early_job_t, 1);
	job->r0 = r0, job->n_rec = n_rec, job->n_frag = e->n_f;
	job->seg_off = kom_calloc(int32_t, 2 * e->n_f + n_rec);
	job->seg_cnt = job->seg_off + e->n_f;
	job->n_hit = job->seg_cnt + e->n_f;
	for (f = 0; f < e->n_f; ++f)
		job->seg_off[f] = e->f_off[f] - r0, job->seg_cnt[f] = e->f_cnt[f];
	job->hit = kom_calloc(mb_hit_t*, n_rec);
	job->seq = kom_malloc(mb_bseq1_t, n_rec);
	for (r = 0; r < n_rec; ++r)
		job->seq[r] = *mb_early_rec(e, r0 + r);
	pthread_mutex_lock(&e->mtx);
	kom_grow(mb_early_job_t*, e->job, e->n_job, e->m_job);
	e->job[e->n_job++] = job;
	pthread_cond_signal(&e->cv);
	pthread_mutex_unlock(&e->mtx);
}

static void mb_early_group(mb_early_t *e, int32_t n) // the main thread: groups records n_proc..n-1 as step 0 does
{
	int32_t i, k;
	for (i = e->n_proc > 1? e->n_proc : 1; i < n; ++i) {
		if (e->frag_mode && mb_qname_same(mb_early_rec(e, i-1)->name, mb_early_rec(e, i)->name)) continue; // as step 0
		if (i - e->j0 > 2) { // step 0 aborts on such a fragment
			e->broken = 1;
			return;
		}
		if (e->n_f == e->m_f) {
			e->m_f += (e->m_f>>1) + 16;
			e->f_off = kom_realloc(int32_t, e->f_off, e->m_f);
			e->f_cnt = kom_realloc(int32_t, e->f_cnt, e->m_f);
		}
		e->f_off[e->n_f] = e->j0, e->f_cnt[e->n_f++] = i - e->j0;
		for (k = e->j0; k < i; ++k)
			e->sb_len += mb_early_rec(e, k)->l_seq;
		++e->n_frag; // fragments 0..n_frag-1 are closed and fragment n_frag starts at record i
		if (e->sb_len >= e->opt->sb_len || e->n_frag - e->sb_off >= e->opt->sb_seq) { // as step 0 at fragment n_frag
			mb_early_push(e);
			e->sb_len = 0, e->sb_off = e->n_frag, e->n_f = 0;
		}
		e->j0 = i;
	}
	e->n_proc = n;
}

static void mb_early_map(mb_early_t *e, mb_early_job_t *job, int32_t tid)
{
	step_t v;
	int32_t f, zero = 0;
	int64_t l_buf = 0;
	char *buf = 0;
	for (f = 0; f < 2; ++f) { // f == 0: measure; f == 1: copy the names step 0 trims ("/[12]"), trimmed
		int32_t i;
		if (f == 1) {
			if (l_buf == 0) break;
			buf = kom_malloc(char, l_buf), l_buf = 0;
		}
		for (i = 0; i < job->n_frag; ++i) {
			mb_bseq1_t *t0 = &job->seq[job->seg_off[i]], *t1 = t0 + 1;
			int32_t l0, l1;
			if (job->seg_cnt[i] != 2) continue;
			l0 = strlen(t0->name), l1 = strlen(t1->name);
			if (l0 < 3 || l0 != l1) continue;
			if (t0->name[l0 - 1] != t1->name[l1 - 1] && t0->name[l0 - 2] == '/') {
				if (f == 1) {
					char *c0 = &buf[l_buf], *c1 = &buf[l_buf + l0 - 1];
					memcpy(c0, t0->name, l0 - 2), c0[l0 - 2] = 0, t0->name = c0;
					memcpy(c1, t1->name, l1 - 2), c1[l1 - 2] = 0, t1->name = c1;
				}
				l_buf += 2 * (l0 - 1);
			}
		}
	}
	memset(&v, 0, sizeof(step_t));
	v.p = e->pl, v.n_seq = job->n_rec, v.n_frag = job->n_frag, v.n_sb = 1;
	v.n_hit = job->n_hit, v.seg_off = job->seg_off, v.seg_cnt = job->seg_cnt;
	v.sb_off = &zero, v.sb_cnt = &job->n_frag;
	v.seq = job->seq, v.hit = job->hit, v.tbuf = e->tbuf;
	worker_for_se_batch(&v, 0, tid);
	free(buf);
	free(job->seq);
	job->seq = 0, job->done = 1;
}

static void *mb_early_worker(void *data)
{
	mb_early_worker_t *w = (mb_early_worker_t*)data;
	mb_early_t *e = w->e;
	for (;;) {
		mb_early_job_t *job = 0;
		pthread_mutex_lock(&e->mtx);
		while (e->reading && e->next == e->n_job)
			pthread_cond_wait(&e->cv, &e->mtx);
		if (e->reading) job = e->job[e->next++]; // once the read ends, step 1 maps what is left
		pthread_mutex_unlock(&e->mtx);
		if (job == 0) break;
		mb_early_map(e, job, w->tid);
	}
	return 0;
}

static void mb_early_end(mb_early_t *e) // the reader is done
{
	pthread_mutex_lock(&e->mtx);
	e->reading = 0;
	pthread_cond_broadcast(&e->cv);
	pthread_cond_broadcast(&e->cv_pub);
	pthread_mutex_unlock(&e->mtx);
}

static void mb_early_run(mb_early_t *e, const pipeline_t *pl) // the index has loaded; returns when the read ends
{
	int32_t i, go;
	e->pl = pl;
	pthread_mutex_lock(&e->mtx);
	if ((go = e->reading) != 0)
		__atomic_store_n(&e->want, 1, __ATOMIC_RELEASE);
	pthread_mutex_unlock(&e->mtx);
	if (!go) return;
	e->n_tbuf = pl->opt->n_thread - 1;
	e->tbuf = kom_calloc(mb_tbuf_t*, e->n_tbuf);
	e->tid = kom_calloc(pthread_t, e->n_tbuf);
	e->w = kom_calloc(mb_early_worker_t, e->n_tbuf);
	for (i = 0; i < e->n_tbuf; ++i)
		e->tbuf[i] = mb_tbuf_init(pl->opt->flag & MB_F_NO_KALLOC);
	for (i = 0; i < e->n_tbuf; ++i) {
		e->w[i].e = e, e->w[i].tid = i;
		if (pthread_create(&e->tid[i], 0, mb_early_worker, &e->w[i]) != 0) break;
	}
	e->n_worker = i;
	for (;;) { // group the records the reader publishes until the read ends
		int32_t n = 0, reading;
		pthread_mutex_lock(&e->mtx);
		while ((reading = e->reading) != 0 && (n = __atomic_load_n(&e->n_pub, __ATOMIC_ACQUIRE)) <= e->n_proc)
			pthread_cond_wait(&e->cv_pub, &e->mtx);
		pthread_mutex_unlock(&e->mtx);
		if (!reading) break;
		mb_early_group(e, n);
		if (e->broken) break;
	}
}

static void mb_early_join(mb_early_t *e)
{
	int32_t i;
	for (i = 0; i < e->n_worker; ++i)
		pthread_join(e->tid[i], 0);
	for (i = 0; i < e->n_tbuf; ++i)
		mb_tbuf_destroy(e->tbuf[i]);
	e->n_worker = e->n_tbuf = 0;
}

static int32_t mb_early_merge(pipeline_t *p, step_t *s) // step 1 of the first batch; returns the number of sub-batches left to map, moved to the front of s->sb_off[] and s->sb_cnt[]
{
	mb_early_t *e = (mb_early_t*)p->early;
	int32_t k, f, r, m, ok = 1;
	p->early = 0;
	for (k = 0; k < e->n_job && ok; ++k) { // the queued sub-batches are step 0's first ones, in order
		const mb_early_job_t *job = e->job[k];
		if (!job->done) continue;
		if (k >= s->n_sb || s->sb_cnt[k] != job->n_frag) ok = 0;
		for (f = 0; ok && f < job->n_frag; ++f)
			if (s->seg_off[s->sb_off[k] + f] != job->r0 + job->seg_off[f] || s->seg_cnt[s->sb_off[k] + f] != job->seg_cnt[f])
				ok = 0;
	}
	for (k = m = 0; k < s->n_sb; ++k)
		if (!ok || k >= e->n_job || !e->job[k]->done)
			s->sb_off[m] = s->sb_off[k], s->sb_cnt[m++] = s->sb_cnt[k];
	for (k = 0; ok && k < e->n_job; ++k) {
		mb_early_job_t *job = e->job[k];
		if (!job->done) continue;
		for (r = 0; r < job->n_rec; ++r)
			s->n_hit[job->r0 + r] = job->n_hit[r], s->hit[job->r0 + r] = job->hit[r];
		job->taken = 1;
	}
	mb_early_destroy(e);
	return m;
}

typedef struct {
	int32_t n_fp, with_qual, with_comment, frag_mode;
	int64_t mb_size, max_mb_size;
	const char **fn;
	mb_bseq_file_t **fp; // 0 when fn[err_i] failed to open, with errno err_no
	int32_t err_i, err_no;
	mb_bseq1_t *seq; // the first batch
	int n_seq;
	kstring_t warn; // bseq's warnings while reading it
	pthread_t tid;
	mb_early_t *early; // 0 when the batch is not mapped while it is read
} mb_prefetch_t;

static void mb_prefetch_read(mb_prefetch_t *q)
{
	int32_t i, j;
	q->fp = kom_calloc(mb_bseq_file_t*, q->n_fp);
	for (i = 0; i < q->n_fp; ++i) { // as mb_open_bseqs()
		if ((q->fp[i] = mb_bseq_open(q->fn[i])) == 0) {
			q->err_no = errno, q->err_i = i;
			for (j = 0; j < i; ++j)
				mb_bseq_close(q->fp[j]);
			free(q->fp);
			q->fp = 0;
			return;
		}
	}
	for (i = 0; i < q->n_fp; ++i) // inflate each file on its own thread while this one parses
		mb_bseq_ahead(q->fp[i], 1);
	if (q->n_fp > 1) q->seq = mb_bseq_read_frag(q->n_fp, q->fp, q->mb_size, q->with_qual, q->with_comment, &q->n_seq);
	else q->seq = mb_bseq_read(q->fp[0], q->mb_size, q->with_qual, q->with_comment, q->frag_mode, 40000, q->max_mb_size, &q->n_seq); // 40000: min_read_cnt in worker_pipeline()
	for (i = 0; i < q->n_fp; ++i) // stop the inflate threads: the thread that parses the later batches inflates them
		mb_bseq_ahead(q->fp[i], 0);
}

static void *mb_prefetch_worker(void *data)
{
	mb_prefetch_t *q = (mb_prefetch_t*)data;
	mb_prefetch_read(q);
	if (q->early) mb_early_end(q->early);
	return 0;
}

static int mb_prefetch_safe(int32_t n, const char **fn, const char *fn_out)
{
	struct stat so, si;
	int32_t i, has_out;
	has_out = fn_out && strcmp(fn_out, "-") != 0 && stat(fn_out, &so) == 0;
	for (i = 0; i < n; ++i) {
		if ((fn[i] && strcmp(fn[i], "-")? stat(fn[i], &si) : fstat(0, &si)) != 0) return 0;
		if (has_out && si.st_dev == so.st_dev && si.st_ino == so.st_ino) return 0;
	}
	return 1;
}

static mb_prefetch_t *mb_prefetch_start(const mb_opt_t *opt, int32_t n, const char **fn, const char *fn_out)
{
	mb_prefetch_t *q;
	if (n < 1 || opt->n_thread <= 1 || !mb_prefetch_safe(n, fn, fn_out)) return 0;
	q = kom_calloc(mb_prefetch_t, 1);
	q->n_fp = n, q->fn = fn;
	q->with_qual = !(opt->flag & MB_F_PAF);
	q->with_comment = !!(opt->flag & MB_F_COPY_COMMENT);
	q->frag_mode = (n > 1 || !!(opt->flag & MB_F_PE));
	q->mb_size = opt->mb_size, q->max_mb_size = opt->max_mb_size;
	mb_bseq_warn_buffer(&q->warn);
	if (kom_dbg_flag == 0 && opt->sb_len > 0 && opt->sb_seq > 0) {
		q->early = mb_early_init(opt, q->frag_mode);
		mb_bseq_read_hook(mb_early_hook, q->early);
	}
	if (pthread_create(&q->tid, 0, mb_prefetch_worker, q) != 0) {
		mb_bseq_warn_buffer(0);
		mb_bseq_read_hook(0, 0);
		if (q->early) mb_early_destroy(q->early);
		free(q);
		return 0;
	}
	return q;
}

static int32_t mb_map_file_pre(const mb_opt_t *opt, const mb_idx_t *idx, int32_t n, const char **fn, const char *hdr, const char *fn_out, mb_prefetch_t *q)
{
	int32_t i, pl_thread;
	pipeline_t pl;
	if (n < 1) return -1;
	memset(&pl, 0, sizeof(pipeline_t));
	pl.fp_out = fn_out == 0 || strcmp(fn_out, "-") == 0? stdout : fopen(fn_out, "wb");
	if (pl.fp_out == 0) return -1;
	if (hdr) fputs(hdr, pl.fp_out);
	pl.n_fp = n;
	if (q) {
		if (q->early) {
			pl.opt = opt, pl.idx = idx;
			mb_early_run(q->early, &pl);
		}
		pthread_join(q->tid, 0);
		mb_bseq_warn_buffer(0);
		mb_bseq_read_hook(0, 0);
		if (q->early) {
			mb_early_join(q->early);
			if (q->early->n_job > 0) pl.early = q->early;
			else mb_early_destroy(q->early);
		}
		if (q->fp == 0 && kom_verbose >= 1)
			fprintf(stderr, "ERROR: failed to open file '%s': %s\n", fn[q->err_i], strerror(q->err_no));
		if (q->fp && q->warn.l > 0)
			fwrite(q->warn.s, 1, q->warn.l, stderr);
		pl.fp = q->fp;
		pl.pre = 1, pl.pre_seq = q->seq, pl.pre_n = q->n_seq;
		free(q->warn.s);
		free(q);
	} else pl.fp = mb_open_bseqs(pl.n_fp, fn);
	if (pl.fp == 0) return -1;
	pl.opt = opt, pl.idx = idx;
	pl.mb_size = opt->mb_size;
	pl_thread = opt->n_thread <= 2? opt->n_thread : 3;

	kt_pipeline(pl_thread, worker_pipeline, &pl, 3);
	if (pl.early) mb_early_destroy((mb_early_t*)pl.early); // the first batch was empty

	if (pl.fp_out != stdout) fclose(pl.fp_out);
	for (i = 0; i < n; ++i)
		mb_bseq_close(pl.fp[i]);
	free(pl.fp);
	return 0;
}

int32_t mb_map_file(const mb_opt_t *opt, const mb_idx_t *idx, int32_t n, const char **fn, const char *hdr, const char *fn_out)
{
	return mb_map_file_pre(opt, idx, n, fn, hdr, fn_out, 0);
}

/********************
 * Arbitrary header *
 ********************/

static void mb_insert_hdr(kstring_t *out, const char *str)
{
	out->l = 0;
	if (str[0] == '@' && strlen(str) > 4) {
		char *dup;
		dup = mb_escape(kom_strdup(str));
		kom_sprintf_lite(out, "%s\n", dup);
		free(dup);
	} else if (str[0] != '@') {
		int dret;
		kstring_t tmp = {0,0,0};
		gzFile fp;
		kstream_t *ks;
		fp = gzopen(str, "r");
		if (fp == 0) return;
		ks = ks_init(fp);
		while (ks_getuntil(ks, KS_SEP_LINE, &tmp, &dret) >= 0)
			if (tmp.s[0] == '@' && tmp.l > 4 && tmp.s[3] == '\t')
				kom_sprintf_lite(out, "%s\n", tmp.s);
		free(tmp.s);
		ks_destroy(ks);
		gzclose(fp);
	}
}

/*******
 * CLI *
 *******/

static ko_longopt_t long_options[] = {
	{ "kalloc",       ko_no_argument,       301 },
	{ "outn",         ko_required_argument, 302 },
	{ "rescue",       ko_required_argument, 304 },
	{ "eqx",          ko_no_argument,       305 },
	{ "pe",           ko_required_argument, 306 },
	{ "long",         ko_optional_argument, 307 },
	{ "adap",         ko_required_argument, 308 },
	{ "chain-only",   ko_no_argument,       309 },
	{ "meth",         ko_no_argument,       310 },
	{ "hic",          ko_no_argument,       311 },
	{ "xa",           ko_required_argument, 312 },
	{ "mmap",         ko_optional_argument, 313 },
	{ "xa-ratio",     ko_required_argument, 314 },
	{ "outs",         ko_required_argument, 315 },
	{ "dbg-aln-seq",  ko_no_argument,       601 },
	{ "dbg-anchor",   ko_no_argument,       602 },
	{ "dbg-seed",     ko_no_argument,       603 },
	{ "dbg-qname",    ko_no_argument,       604 },
	{ "dbg-aln-pe",   ko_no_argument,       605 },
	{ "dbg-an-pos",   ko_no_argument,       606 }, // anchor position
	{ "version",      ko_no_argument,       901 },
	{ "help",         ko_no_argument,       902 },
	{ 0, 0, 0 }
};

static int usage_map(FILE *fp, const mb_opt_t *opt)
{
	fprintf(fp, "Usage: minibwa map [options] <in.idx> <in.fastq>\n");
	fprintf(fp, "Options:\n");
	fprintf(fp, "  Common:\n");
	fprintf(fp, "    -f               output PAF (SAM by default)\n");
	fprintf(fp, "    -t INT           number of worker threads [%d]\n", opt->n_thread);
	fprintf(fp, "    -l NUM           treat reads <NUM as short reads in the default adaptive mode [%d]\n", opt->max_sr_len);
	fprintf(fp, "    -R STR           SAM read group line in a format like '@RG\\tID:foo\\tSM:bar' []\n");
	fprintf(fp, "    -b STR           output a base alignment tag: cs, ds or MD []\n");
	fprintf(fp, "    --hic            map Hi-C reads; equivalent to option -5P\n");
	fprintf(fp, "    --meth           map *directional* bisulfite sequencing reads\n");
	fprintf(fp, "  Mapping:\n");
	fprintf(fp, "    -k INT           min seed length [%d]\n", opt->min_len);
	fprintf(fp, "    -c NUM           max seed occurrences [%d]\n", opt->max_occ);
	fprintf(fp, "    -g NUM           max gap size, controlling extension and chain breaking [%d]\n", opt->max_gap);
	fprintf(fp, "    -w NUM           bandwidth [%d]\n", opt->bw);
	fprintf(fp, "    -W NUM           long bandwidth (for long reads or the adaptive mode) [%d]\n", opt->bw_long);
    fprintf(fp, "    -m INT           min chaining score [%d]\n", opt->min_chain_score);
	fprintf(fp, "    -p FLOAT         min secondary-to-primary score ratio [%g]\n", opt->pri_ratio);
	fprintf(fp, "    -N INT           retain at most INT secondary alignments [%d]\n", opt->best_n);
	fprintf(fp, "    --chain-only     perform chaining only without base alignment; force -fP\n");
	fprintf(fp, "    -x STR           preset (sr, lr or adap for mixed short/long reads) [adap]\n");
	fprintf(fp, "  Alignment:\n");
	fprintf(fp, "    -A INT           matching score [%d]\n", opt->a);
	fprintf(fp, "    -B INT           mismatching openalty [%d]\n", opt->b);
	fprintf(fp, "    -O INT1[,INT2]   gap open penalty [%d,%d]\n", opt->q, opt->q2);
	fprintf(fp, "    -E INT1[,INT2]   gap extension penalty [%d,%d]\n", opt->e, opt->e2);
	fprintf(fp, "    -s INT           suppress alignment with DP score lower than INT*{-A} [%d]\n", opt->min_dp_max);
	fprintf(fp, "  Paired-end:\n");
	fprintf(fp, "    -P               skip pairing and mate rescue\n");
	fprintf(fp, "    --rescue=INT     mate rescue for up to INT candidates; 0 to skip rescue [%d]\n", opt->max_rescue);
	fprintf(fp, "    -I INT[,INT[,INT[,INT]]]\n");
	fprintf(fp, "                     mean, stddev, max and min of isize distribution [inferred]\n");
	fprintf(fp, "  Input/Output:\n");
	fprintf(fp, "    -o FILE          output file name [stdout]\n");
	fprintf(fp, "    -u               don't output unmapped reads\n");
	fprintf(fp, "    --outn=NUM       output up to {NUM,-N} secondary alignments [0]\n");
	fprintf(fp, "    --outs=FLOAT     output a secondary hit if score at least FLOAT*bestScore [%g]\n", opt->out_s);
	fprintf(fp, "    --xa=NUM         if <=NUM hits with score >%g%% of the best hit, output them to XA [%d]\n", opt->out_s*100.0, opt->xa_max);
	fprintf(fp, "    -y               copy FASTA/Q comments to output\n");
	fprintf(fp, "    -Y               use soft clipping for supplementary alignments\n");
	fprintf(fp, "    -H STR           if STR starts with @, insert to header; or insert lines in file STR []\n");
	fprintf(fp, "    -5               take the alignment with the smallest query position as primary\n");
	fprintf(fp, "    -K NUM1[,NUM2]   process NUM1-NUM2 bp of query sequences in a batch [100m,1g]\n");
	fprintf(fp, "    --mmap[=lite]    load the index via memory mapped files (slower mapping) []\n");
	fprintf(fp, "    --version        print version number\n");
	fprintf(fp, "    --help           print this help message\n");
	return fp == stdout? 0 : 1;
}

static inline void yes_or_no(mb_opt_t *opt, uint64_t flag, int long_idx, const char *arg, int yes_to_set)
{
	if (yes_to_set) {
		if (strcmp(arg, "yes") == 0 || strcmp(arg, "y") == 0) opt->flag |= flag;
		else if (strcmp(arg, "no") == 0 || strcmp(arg, "n") == 0) opt->flag &= ~flag;
		else fprintf(stderr, "[WARNING]\033[1;31m option '--%s' only accepts 'yes' or 'no'.\033[0m\n", long_options[long_idx].name);
	} else {
		if (strcmp(arg, "yes") == 0 || strcmp(arg, "y") == 0) opt->flag &= ~flag;
		else if (strcmp(arg, "no") == 0 || strcmp(arg, "n") == 0) opt->flag |= flag;
		else fprintf(stderr, "[WARNING]\033[1;31m option '--%s' only accepts 'yes' or 'no'.\033[0m\n", long_options[long_idx].name);
	}
}

static void set_ins_size(mb_opt_t *opt, const char *arg)
{
	char *q;
	opt->pe_avg = (int32_t)(.499 + strtod(arg, &q));
	opt->pe_std = (int32_t)(.499 + (*q == ','? strtod(q + 1, &q) : opt->pe_avg * 0.1));
	opt->pe_hi = *q == ','? strtol(q + 1, &q, 10) : opt->pe_avg + opt->pe_std * 4;
	opt->pe_lo = *q == ','? strtol(q + 1, &q, 10) : opt->pe_avg - opt->pe_std * 4;
	if (opt->pe_lo < 1) opt->pe_lo = 1;
	opt->flag |= MB_F_PE_PREDEF;
}

int main_map(int argc, char *argv[])
{
	const char *opt_str = "x:o:k:c:m:p:A:B:U:b:O:E:t:K:N:PyYR:H:aul:w:W:g:5s:fI:";
	int32_t c, use_mmap = 0, mmap_preload = 1, is_meth = 0;
	mb_idx_t *idx;
	mb_prefetch_t *pre;
	mb_opt_t mo;
	char *fn_out = 0, *rg_line = 0, *s;
	ketopt_t o = KETOPT_INIT;
	kstring_t hdr_ins = {0,0,0}, hdr = {0,0,0};

	mb_opt_init(&mo);
	while ((c = ketopt(&o, argc, argv, 1, opt_str, long_options)) >= 0) { // test command line options and apply option -x/preset first
		if (c == 'x') {
			if (mb_opt_preset(&mo, o.arg) < 0) {
				fprintf(stderr, "[ERROR] unknown preset '%s'\n", o.arg);
				return 1;
			}
		} else if (c == ':') {
			fprintf(stderr, "[ERROR] missing option argument\n");
			return 1;
		} else if (c == '?') {
			fprintf(stderr, "[ERROR] unknown option in \"%s\"\n", argv[o.erri]);
			return 1;
		}
	}
	o = KETOPT_INIT;
	while ((c = ketopt(&o, argc, argv, 1, opt_str, long_options)) >= 0) {
		if (c == 'k') mo.min_len = atoi(o.arg);
		else if (c == 'c') mo.max_occ = kom_parse_num(o.arg, 0);
		else if (c == 'p') mo.pri_ratio = atof(o.arg);
		else if (c == 'm') mo.min_chain_score = atoi(o.arg);
		else if (c == 'N') mo.best_n = atoi(o.arg);
		else if (c == 'A') mo.a = atoi(o.arg);
		else if (c == 'B') mo.b = atoi(o.arg);
		else if (c == 'U') mo.pen_unpair = atoi(o.arg);
		else if (c == 'l') mo.max_sr_len = kom_parse_num(o.arg, 0);
		else if (c == 'g') mo.max_gap = kom_parse_num(o.arg, 0);
		else if (c == 'w') mo.bw = kom_parse_num(o.arg, 0);
		else if (c == 'W') mo.bw_long = kom_parse_num(o.arg, 0);
		else if (c == 'a') mo.flag &= ~MB_F_PAF;
		else if (c == 'f') mo.flag |= MB_F_PAF;
		else if (c == 'u') mo.flag |= MB_F_NO_UNMAP;
		else if (c == 'y') mo.flag |= MB_F_COPY_COMMENT;
		else if (c == 'Y') mo.flag |= MB_F_SUPP_SOFT;
		else if (c == '5') mo.flag |= MB_F_PRIMARY5;
		else if (c == 'P') mo.flag |= MB_F_NO_PAIRING;
		else if (c == 's') mo.min_dp_max = atoi(o.arg);
		else if (c == 'o') fn_out = o.arg;
		else if (c == 't') mo.n_thread = atoi(o.arg);
		else if (c == 'R') rg_line = o.arg;
		else if (c == 'H') mb_insert_hdr(&hdr_ins, o.arg);
		else if (c == 'I') set_ins_size(&mo, o.arg);
		else if (c == 301) { // --kalloc
			yes_or_no(&mo, MB_F_NO_KALLOC, o.longidx, o.arg, 0);
		} else if (c == 302) { // --outn
			mo.out_n = kom_parse_num(o.arg, 0);
		} else if (c == 304) { // --rescue
			mo.max_rescue = atoi(o.arg);
		} else if (c == 305) { // --eqx
			mo.flag |= MB_F_EQX;
		} else if (c == 306) { // --pe
			yes_or_no(&mo, MB_F_PE, o.longidx, o.arg, 1);
		} else if (c == 307) { // --long
			if (o.arg == 0) mo.flag |= MB_F_LONG;
			else yes_or_no(&mo, MB_F_LONG, o.longidx, o.arg, 1);
		} else if (c == 308) { // --adap
			yes_or_no(&mo, MB_F_ADAP, o.longidx, o.arg, 1);
		} else if (c == 309) { // --chain-only
			mo.flag |= MB_F_NO_ALN;
		} else if (c == 310) { // --meth
			mo.flag |= MB_F_METH;
		} else if (c == 311) { // --hic
			mo.flag |= MB_F_PRIMARY5 | MB_F_NO_PAIRING;
		} else if (c == 312) { // --xa
			mo.xa_max = kom_parse_num(o.arg, 0);
		} else if (c == 313) { // --mmap
			use_mmap = 1;
			if (o.arg != 0 && strcmp(o.arg, "lite") == 0) mmap_preload = 0;
		} else if (c == 314 || c == 315) { // --outs or --xa-ratio
			mo.out_s = atof(o.arg);
		} else if (c == 601) { // --dbg-aln-seq
			kom_dbg_flag |= MB_DBG_ALN_SEQ;
		} else if (c == 602) { // --dbg-anchor
			kom_dbg_flag |= MB_DBG_ANCHOR;
		} else if (c == 603) { // --dbg-seed
			kom_dbg_flag |= MB_DBG_SEED;
		} else if (c == 604) { // --dbg-qname
			kom_dbg_flag |= MB_DBG_QNAME;
		} else if (c == 605) { // --dbg-aln-pe
			kom_dbg_flag |= MB_DBG_ALN_PE;
		} else if (c == 606) { // --dbg-an-pos
			kom_dbg_flag |= MB_DBG_AN_POS;
		} else if (c == 'K') {
			mo.mb_size = mo.max_mb_size = kom_parse_num(o.arg, &s);
			if (*s == ',') mo.max_mb_size = kom_parse_num(s + 1, &s);
			if (mo.max_mb_size < mo.mb_size)
				mo.max_mb_size = mo.mb_size;
		} else if (c == 'O') {
			mo.q = mo.q2 = strtol(o.arg, &s, 10);
			if (*s == ',') mo.q2 = strtol(s + 1, &s, 10);
		} else if (c == 'E') {
			mo.e = mo.e2 = strtol(o.arg, &s, 10);
			if (*s == ',') mo.e2 = strtol(s + 1, &s, 10);
		} else if (c == 'b') {
			mo.flag &= ~(MB_F_WRITE_CS|MB_F_WRITE_DS|MB_F_WRITE_MD);
			if (strcmp(o.arg, "cs") == 0) mo.flag |= MB_F_WRITE_CS;
			else if (strcmp(o.arg, "ds") == 0) mo.flag |= MB_F_WRITE_DS;
			else if (strcmp(o.arg, "MD") == 0 || strcmp(o.arg, "md") == 0) mo.flag |= MB_F_WRITE_MD;
			else if (kom_verbose >= 2) {
				mo.flag |= MB_F_WRITE_CS;
				fprintf(stderr, "[WARNING]\033[1;31m -b only takes 'cs', 'ds' or 'MD'. Invalid values are assumed to be 'cs'.\033[0m\n");
			}
		} else if (c == 901) { // --version
			puts(MB_VERSION);
			exit(0);
		} else if (c == 902) { // --help
			return usage_map(stdout, &mo);
		}
	}
	if (mo.flag & MB_F_NO_ALN) mo.flag |= MB_F_NO_PAIRING | MB_F_PAF;
	if (argc - o.ind < 2)
		return usage_map(stderr, &mo);

	is_meth = !!(mo.flag & MB_F_METH);
	pre = mb_prefetch_start(&mo, argc - (o.ind + 1), (const char**)&argv[o.ind+1], fn_out);
	idx = use_mmap? mb_idx_load_mmap(argv[o.ind], is_meth, mmap_preload) : mb_idx_load_mt(argv[o.ind], is_meth, mo.n_thread, kt_for);
	kom_assert(idx, "failed to load the index.");
	if (kom_verbose >= 3)
		fprintf(stderr, "[M::%s::%.3f*%.2f] index loaded\n", __func__, kom_realtime(), kom_percent_cpu());

	if (!(mo.flag & MB_F_PAF)) {
		int ret;
		ret = mb_fmt_sam_hdr(&hdr, idx->l2b, rg_line, MB_VERSION, argc, argv);
		if (ret < 0) return 1; // TODO: free idx and out.s
		if (hdr_ins.l > 0) kom_sprintf_lite(&hdr, "%s", hdr_ins.s);
	}
	if (hdr_ins.s) free(hdr_ins.s);

	mb_map_file_pre(&mo, idx, argc - (o.ind + 1), (const char**)&argv[o.ind+1], hdr.s, fn_out, pre);
	free(hdr.s);
	mb_idx_destroy(idx);
	return 0;
}

static int usage_mem(FILE *fp, const mb_opt_t *opt)
{
	fprintf(fp, "Usage: minibwa mem [options] <in.idx> <in.fastq>\n");
	fprintf(fp, "Options:\n");
	fprintf(fp, "  Algorithm:\n");
	fprintf(fp, "    -t INT         number of threads [%d]\n", opt->n_thread);
	fprintf(fp, "    -k INT         minimum seed length [%d]\n", opt->min_len);
	fprintf(fp, "    -w NUM         bandwidth [%d]\n", opt->bw);
	fprintf(fp, "    *r FLOAT       inside seed ratio []\n");
	fprintf(fp, "    *y INT         seed occurrence for the 3rd round seeding []\n");
	fprintf(fp, "    -c NUM         max seed occurrences [%d]\n", opt->max_occ);
	fprintf(fp, "    -D FLOAT       mask level [%g]\n", opt->mask_level);
	fprintf(fp, "    -W INT         min chaining score [%d]\n", opt->min_chain_score);
	fprintf(fp, "    -m INT         mate rescue for up to INT candidates; 0 to skip rescue [%d]\n", opt->max_rescue);
	fprintf(fp, "    -S             skip mate rescue\n");
	fprintf(fp, "    -P             skip pairing\n");
	fprintf(fp, "  Scoring:\n");
	fprintf(fp, "    *A INT         matching score []\n");
	fprintf(fp, "    *B INT         mismatching penalty []\n");
	fprintf(fp, "    *O INT         gap open penalty []\n");
	fprintf(fp, "    *E INT         gap extension penalty []\n");
	fprintf(fp, "    *L INT         end bonus []\n");
	fprintf(fp, "    *U INT         penalty for an unpaired read pair []\n");
	fprintf(fp, "    *d INT         off-diagonal X-dropoff []\n");
	fprintf(fp, "    *x STR         preset []\n");
	fprintf(fp, "  Input/output:\n");
	fprintf(fp, "    *p             smart pairing\n");
	fprintf(fp, "    -R STR         SAM read group line in a format like '@RG\\tID:foo\\tSM:bar' []\n");
	fprintf(fp, "    -H STR         if STR starts with @, insert to header; or insert lines in file STR []\n");
	fprintf(fp, "    -o FILE        output file name [stdout]\n");
	fprintf(fp, "    *j             treat ALT contigs as part of the primary assembly\n");
	fprintf(fp, "    -5             take the alignment with the smallest query position as primary\n");
	fprintf(fp, "    *q             don't modify mapQ of supplementary alignments\n");
	fprintf(fp, "    *K NUM         batch size []\n");
	fprintf(fp, "    -v INT         verbose level: 1=error, 2=warning, 3=message, 4+=debugging [%d]\n", kom_verbose);
	fprintf(fp, "    -T INT         suppress alignment with DP score lower than INT*{-A} [%d]\n", opt->min_dp_max);
	fprintf(fp, "    -h INT         max secondary alignments at the XA tag [%d]\n", opt->xa_max);
	fprintf(fp, "    -z FLOAT       skip XA if score < FLOAT*bestScore [%g]\n", opt->out_s);
	fprintf(fp, "    -a             output all alignments\n");
	fprintf(fp, "    -C             copy FASTA/Q comments to output\n");
	fprintf(fp, "    *V             output the reference FASTA header in the XR tag\n");
	fprintf(fp, "    -Y             use soft clipping for supplementary alignments\n");
	fprintf(fp, "    *M             mark shorter split hits as secondary\n");
	fprintf(fp, "    -I INT[,INT[,INT[,INT]]]\n");
	fprintf(fp, "                   mean, stddev, max and min of isize distribution [inferred]\n");
	fprintf(fp, "    *u             output XB instead of XA\n");
	fprintf(fp, "Notes:\n");
	fprintf(fp, "  - \"minibwa mem\" aims to match the \"bwa mem\" command-line interface\n");
	fprintf(fp, "  - '*' options are ignored as they are missing or incompatible with minibwa\n");
	fprintf(fp, "  - minibwa and bwa-mem may output different alignments\n");
	return fp == stdout? 0 : 1;
}

int main_mem(int argc, char *argv[])
{
	int32_t c, ret;
	ketopt_t o = KETOPT_INIT;
	mb_opt_t mo;
	char *fn_out = 0, *rg_line = 0;
	kstring_t hdr_ins = {0,0,0}, hdr = {0,0,0};
	mb_idx_t *idx;
	mb_prefetch_t *pre;

	mb_opt_init(&mo);
	mo.flag |= MB_F_WRITE_MD; // bwa-mem always writes MD
	while ((c = ketopt(&o, argc, argv, 1, "t:k:w:d:r:y:c:D:W:m:SPA:B:O:E:L:U:x:pR:H:o:j5qK:v:T:h:z:aCVYMuI:", 0)) >= 0) {
		// algorithm
		if (c == 't') mo.n_thread = atoi(o.arg);
		else if (c == 'k') mo.min_len = atoi(o.arg);
		else if (c == 'w') mo.bw = kom_parse_num(o.arg, 0);
		else if (c == 'r' || c == 'y') {}
		else if (c == 'W') mo.min_chain_score = atoi(o.arg);
		else if (c == 'c') mo.max_occ = kom_parse_num(o.arg, 0);
		else if (c == 'D') mo.mask_level = atof(o.arg);
		else if (c == 'm') mo.max_rescue = atoi(o.arg);
		else if (c == 'S') mo.max_rescue = 0;
		else if (c == 'P') mo.flag |= MB_F_NO_PAIRING;
		// scoring
		else if (c == 'A' || c == 'B' || c == 'O' || c == 'E' || c == 'L' || c == 'U' || c == 'd' || c == 'x') {}
		// input/output
		else if (c == 'p' || c == 'j' || c == 'q' || c == 'K' || c == 'V' || c == 'M' || c == 'u') {}
		else if (c == 'a') mo.out_n = 1000000;
		else if (c == 'h') mo.xa_max = atoi(o.arg);
		else if (c == 'z') mo.out_s = atof(o.arg);
		else if (c == 'R') rg_line = o.arg;
		else if (c == 'H') mb_insert_hdr(&hdr_ins, o.arg);
		else if (c == 'o') fn_out = o.arg;
		else if (c == '5') mo.flag |= MB_F_PRIMARY5;
		else if (c == 'v') kom_verbose = atoi(o.arg);
		else if (c == 'T') mo.min_dp_max = atoi(o.arg);
		else if (c == 'C') mo.flag |= MB_F_COPY_COMMENT;
		else if (c == 'Y') mo.flag |= MB_F_SUPP_SOFT;
		else if (c == 'I') set_ins_size(&mo, o.arg);
	}
	if (argc - o.ind < 2) return usage_mem(stderr, &mo);

	pre = mb_prefetch_start(&mo, argc - (o.ind + 1), (const char**)&argv[o.ind+1], fn_out);
	idx = mb_idx_load_mt(argv[o.ind], !!(mo.flag & MB_F_METH), mo.n_thread, kt_for);
	kom_assert(idx, "failed to load the index.");
	if (kom_verbose >= 3)
		fprintf(stderr, "[M::%s::%.3f*%.2f] index loaded\n", __func__, kom_realtime(), kom_percent_cpu());

	ret = mb_fmt_sam_hdr(&hdr, idx->l2b, rg_line, MB_VERSION, argc, argv);
	if (ret < 0) return 1; // TODO: free idx and out.s
	if (hdr_ins.l > 0) kom_sprintf_lite(&hdr, "%s", hdr_ins.s);
	if (hdr_ins.s) free(hdr_ins.s);
	mb_map_file_pre(&mo, idx, argc - (o.ind + 1), (const char**)&argv[o.ind+1], hdr.s, fn_out, pre);
	free(hdr.s);
	mb_idx_destroy(idx);
	return 0;
}
