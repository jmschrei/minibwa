#include <string.h>
#include <assert.h>
#include <stdio.h>
#include "mbpriv.h"
#include "kalloc.h"
#include "ksort.h"

#define key_sai0(a) ((a).x[0])
KRADIX_SORT_INIT(mb_sai0, mb_sai_t, key_sai0, 8)

#define key_sais(a) ((a).size)
KRADIX_SORT_INIT(mb_sais, mb_sai_t, key_sais, 8)

#define key_saii(a) ((a).info)
KRADIX_SORT_INIT(mb_saii, mb_sai_t, key_saii, 8)

#define key_anchor(a) ((a).tpos)
KRADIX_SORT_INIT(mb_anchor, mb_anchor_t, key_anchor, 8)

/***********
 * Seeding *
 ***********/

void mb_seed_intv(void *km, const mb_bwt_t *bwt, int32_t len, const uint8_t *seq, int32_t min_len, int32_t max_sub_occ, mb_sai_v *v)
{
	int64_t x = 0, i, n_a0;
	mb_sai_t p;

	v->n = 0;
	do { // pass 1: standard SMEMs
		x = mb_bwt_smem(bwt, len, seq, x, min_len, 1, &p);
		if (p.size > 0) {
			Kgrow(km, mb_sai_t, v->a, v->n, v->m);
			v->a[v->n++] = p;
		}
	} while (x < len);

	n_a0 = v->n;
	for (i = 0; i < n_a0; ++i) { // pass 2: sub-SMEMs
		int32_t sub_min_len;
		uint32_t st = v->a[i].info>>32, en = (uint32_t)v->a[i].info;
		if (en - st < min_len * 2 || v->a[i].size > max_sub_occ)
			continue;
		x = st;
		sub_min_len = (en - st) / 2 > min_len? (en - st) / 2 : min_len;
		do { // if two SMEMs have large overlaps, we may find the same sub intervals in both. A rare case not worth optimizing
			x = mb_bwt_smem(bwt, en, seq, x, sub_min_len, v->a[i].size + 1, &p);
			if (p.size > v->a[i].size) {
				Kgrow(km, mb_sai_t, v->a, v->n, v->m);
				v->a[v->n++] = p;
			}
		} while (x < en);
	}
}

typedef struct { // a pending sub-SMEM job
	int32_t min_len, min_occ;
	uint32_t st, en;
	int32_t sid;
} seed_job_t;

typedef struct { // the job source of seed_intv_batch()
	void *km;
	int32_t n_seq, next_seq, min_len, max_sub_occ;
	const int32_t *len;
	uint8_t *const* seq;
	mb_sai_v *v;
	const mb_smem_entry_t *slot;
	int32_t *sid; // read of the job in each slot
	uint8_t *first; // whether the job in each slot is a first pass
	size_t n_job, m_job, i_job;
	seed_job_t *job;
} seed_feed_t;

static int seed_feed(void *data, mb_smem_entry_t *t, int done)
{
	seed_feed_t *f = (seed_feed_t*)data;
	int32_t k = t - f->slot;
	if (done && f->first[k]) { // the first pass of a read has finished; queue its sub-SMEM jobs
		int32_t i = f->sid[k];
		size_t j;
		for (j = 0; j < f->v[i].n; ++j) {
			seed_job_t *p;
			uint32_t st = f->v[i].a[j].info>>32, en = (uint32_t)f->v[i].a[j].info;
			if (en - st < f->min_len * 2 || f->v[i].a[j].size > f->max_sub_occ)
				continue;
			Kgrow(f->km, seed_job_t, f->job, f->n_job, f->m_job);
			p = &f->job[f->n_job++];
			p->min_len = (en - st) / 2 > f->min_len? (en - st) / 2 : f->min_len;
			p->min_occ = f->v[i].a[j].size + 1;
			p->st = st, p->en = en;
			p->sid = i;
		}
	}
	if (f->next_seq < f->n_seq) { // first pass: standard SMEMs
		int32_t i = f->next_seq++;
		t->min_len = f->min_len;
		t->min_occ = 1;
		t->st = 0, t->en = f->len[i];
		t->q = f->seq[i];
		t->v = &f->v[i];
		f->sid[k] = i, f->first[k] = 1;
		return 1;
	}
	if (f->i_job < f->n_job) { // second pass: sub-SMEMs
		const seed_job_t *p = &f->job[f->i_job++];
		t->min_len = p->min_len;
		t->min_occ = p->min_occ;
		t->st = p->st, t->en = p->en;
		t->q = f->seq[p->sid];
		t->v = &f->v[p->sid];
		f->sid[k] = p->sid, f->first[k] = 0;
		return 1;
	}
	return 0;
}

/* mb_seed_intv_batch() with pac, the BWT text's forward strand or NULL, passed
 * to mb_bwt_smem_feed(), which uses it in the first pass (min_occ 1): the x[1]
 * of a size-1 SMEM may then hold the text position of x[0] | MB_SA_KNOWN */
static void seed_intv_batch(void *km, const mb_bwt_t *bwt, const uint64_t *pac, int32_t n_seq, const int32_t *len, uint8_t *const* seq, int32_t min_len, int32_t max_sub_occ, mb_sai_v *v)
{ // identical to mb_seed_intv() though the order of intervals is often different
	const int max_batch_size = 50;
	mb_smem_entry_t *s;
	seed_feed_t f;
	int32_t i;

	// each slot takes the next pending job when its job finishes: the first passes of the
	// reads in order, then the sub-SMEM jobs of reads whose first pass has finished
	for (i = 0; i < n_seq; ++i) v[i].n = 0;
	memset(&f, 0, sizeof(f));
	f.km = km;
	f.n_seq = n_seq, f.min_len = min_len, f.max_sub_occ = max_sub_occ;
	f.len = len, f.seq = seq, f.v = v;
	s = Kcalloc(km, mb_smem_entry_t, max_batch_size);
	f.slot = s;
	f.sid = Kcalloc(km, int32_t, max_batch_size);
	f.first = Kcalloc(km, uint8_t, max_batch_size);
	mb_bwt_smem_feed(km, bwt, pac, max_batch_size, s, seed_feed, &f);
	kfree(km, f.job);
	kfree(km, f.first);
	kfree(km, f.sid);
	kfree(km, s);
}

void mb_seed_intv_batch(void *km, const mb_bwt_t *bwt, int32_t n_seq, const int32_t *len, uint8_t *const* seq, int32_t min_len, int32_t max_sub_occ, mb_sai_v *v)
{
	seed_intv_batch(km, bwt, 0, n_seq, len, seq, min_len, max_sub_occ, v);
}

/*****************************
 * Seed/anchor deduplication *
 *****************************/

static void mb_seed_sort_dedup(mb_sai_v *u)
{
	int64_t i, i0, j;
	if (u->n <= 1) return;
	// sort by ::{x[0],size,info}
	radix_sort_mb_sai0(u->a, u->a + u->n); // sort by ::x[0]
	for (i = 1, i0 = 0; i <= u->n; ++i) {
		if (i == u->n || u->a[i].x[0] != u->a[i0].x[0]) {
			if (i - i0 > 1) {
				int64_t k, k0, n1 = i - i0;
				mb_sai_t *a1 = &u->a[i0];
				radix_sort_mb_sais(&u->a[i0], &u->a[i]); // sort by ::size
				kom_reverse(mb_sai_t, n1, a1);
				for (k = i0 + 1, k0 = i0; k <= i; ++k) {
					if (k == i || u->a[k0].size != u->a[k].size) {
						if (k - k0 > 1)
							radix_sort_mb_saii(&u->a[k0], &u->a[k]); // sort by ::info
						k0 = k;
					}
				}
			}
			i0 = i;
		}
	}
	// dedup
	for (i = 1, j = 0; i < u->n; ++i)
		if (!(u->a[i].x[0] == u->a[j].x[0] && u->a[i].size == u->a[j].size && u->a[i].info == u->a[j].info))
			u->a[++j] = u->a[i];
	u->n = j + 1;
}

/* Remove duplicated anchors. The two-round seeding algorithm may lead an
 * anchor precisely contained in a longer anchor. This routine filters out the
 * shorter anchor. This wouldn't happen to minimap2. */
static void mb_anchor_dedup(mb_anchor_v *v) // NB: assuming sorted by tpos
{
	const int max_back = 100; // to avoid quadratic behavior in the worst case
	int64_t i, j, k;
	for (i = 1; i < v->n; ++i) {
		mb_anchor_t *ai = &v->a[i];
		int64_t tsj, tsi = ai->tpos + 1 - ai->len;
		int32_t qsj, qsi = ai->qpos + 1 - ai->len;
		for (j = i - 1, k = 0; j >= 0 && k < max_back; --j, ++k) {
			mb_anchor_t *aj = &v->a[j];
			if (aj->sid != ai->sid) break;
			if (aj->tpos < tsi) break;
			tsj = aj->tpos + 1 - aj->len;
			qsj = aj->qpos + 1 - aj->len;
			if (tsj >= tsi) { // then j is contained in i
				if (tsj - tsi == qsj - qsi) aj->flt = 1;
			} else if (ai->tpos == aj->tpos) { // then i is contained in j
				if (ai->qpos == aj->qpos) ai->flt = 1;
			}
		}
	}
	for (i = j = 0; i < v->n; ++i)
		if (!v->a[i].flt) v->a[j++] = v->a[i];
	v->n = j;
}

/************************
 * Get contig positions *
 ************************/

#define MB_SA_CHUNK 4096 // max SA rows per mb_bwt_sa_batch() call

/* SA rows mb_anchor() resolves for the SA interval of q: all of them, or up
 * to max_occ sampled ones. Writes them to rows[] unless rows is NULL, and
 * returns their number. */
static int64_t anchor_rows1(const mb_sai_t *q, int32_t max_occ, uint64_t *rows)
{
	int64_t j, m = 0;
	if (q->size <= max_occ) { // get SA for all of them
		if (rows)
			for (j = 0; j < q->size; ++j)
				rows[j] = q->x[0] + j;
		m = q->size;
	} else { // sample up to max_occ
		int32_t n = 0;
		for (j = 0; j < q->size && n < max_occ; ++n) {
			int32_t step = (q->size - j) / (max_occ - n);
			if (step < 1) step = 1;
			if (rows) rows[m] = q->x[0] + j;
			++m;
			j += step;
		}
	}
	return m;
}

/* Sort and deduplicate the seed intervals in u as mb_anchor() does, then
 * append the SA rows it resolves for them, in its order, to (*rows)[*n_rows..].
 * The row of a size-1 interval whose x[1] holds its text position
 * (MB_SA_KNOWN, from mb_bwt_smem_batch_pac()) is replaced by that position,
 * which mb_bwt_sa_batch() passes through.
 *
 * m counts rows in groups, starting a new one when an interval's size would
 * take m past batch_size, and the assert caps a group at m_a rows: it holds
 * for every interval of a valid index and aborts, before any row is resolved,
 * on the corrupt (negative) sizes a single-strand index can give. */
static void anchor_rows(void *km, mb_sai_v *u, int32_t max_occ, int64_t *n_rows, int64_t *m_rows, uint64_t **rows)
{
	const int batch_size = 20;
	int32_t m = 0, m_a = max_occ > batch_size? max_occ : batch_size;
	int64_t i, i0, j, n;
	if (u->n == 0) return;
	mb_seed_sort_dedup(u);
	for (i = 1, i0 = 0; i <= u->n; ++i) {
		if (i == u->n || u->a[i].x[0] != u->a[i0].x[0] || u->a[i].size != u->a[i0].size) {
			n = anchor_rows1(&u->a[i0], max_occ, 0);
			if (u->a[i0].size + m > batch_size) m = 0;
			m += n;
			assert(m <= m_a); // shouldn't happen!
			if (*n_rows + n > *m_rows) {
				*m_rows = *n_rows + n;
				*m_rows += (*m_rows>>1) + 16;
				*rows = Krealloc(km, uint64_t, *rows, *m_rows);
			}
			anchor_rows1(&u->a[i0], max_occ, *rows + *n_rows);
			if (n == 1 && u->a[i0].size == 1)
				for (j = i0; j < i; ++j)
					if (u->a[j].x[1] & MB_SA_KNOWN) {
						(*rows)[*n_rows] = u->a[j].x[1];
						break;
					}
			*n_rows += n;
			i0 = i;
		}
	}
}

// replace the n SA rows in x[] with their SA values
static void anchor_resolve(void *km, const mb_bwt_t *bwt, int64_t n, uint64_t *x)
{
	int64_t i;
	for (i = 0; i < n; i += MB_SA_CHUNK)
		mb_bwt_sa_batch(km, bwt, n - i < MB_SA_CHUNK? n - i : MB_SA_CHUNK, x + i);
}

// add the anchors of seed intervals u->a[st..en), which share the SA interval, at one SA value
static inline void anchor_add(void *km, const mb_idx_t *idx, const mb_sai_v *u, int64_t st, int64_t en, uint64_t sa, int32_t qlen, l2b_meth_t mt, mb_anchor_v *v)
{
	int64_t j;
	for (j = st; j < en; ++j) {
		int32_t qs = u->a[j].info>>32, qe = (int32_t)u->a[j].info;
		int32_t rev, len = qe - qs;
		int64_t tid, cst;
		l2b_meth_t mt_anchor;
		const l2b_ctg_t *ctg;
		mb_anchor_t *q;
		if (mt != L2B_METH_NONE) {
			tid = l2b_intv2cid_meth(idx->l2b, sa, sa + len, &mt_anchor, &cst, &rev);
			if (tid < 0) continue;
			// R1(C2T): keep c2t_f(copy0) and g2a_r(copy2); R2(G2A): keep g2a_f(copy1) and c2t_r(copy3)
			if ((mt_anchor == mt) != (rev == 0)) continue; // filter
		} else {
			tid = l2b_intv2cid(idx->l2b, sa, sa + len, &cst, &rev);
			if (tid < 0) continue;
		}
		rev = !!rev;
		ctg = &idx->l2b->ctg[tid];
		Kgrow(km, mb_anchor_t, v->a, v->n, v->m);
		q = &v->a[v->n++];
		memset(q, 0, sizeof(*q));
		q->sid = tid << 1 | rev;
		q->len = len;
		q->qpos = rev? qlen - 1 - qs : qs + len - 1;
		q->tpos = ctg->off * 2 + ctg->len * rev + cst + len - 1; // for sorting; will be adjusted later
	}
}

static void mb_anchor_split_meth(void *km, const l2b_t *l2b, int32_t min_len, int32_t qlen, const uint8_t *qseq0, l2b_meth_t mt0, mb_anchor_v *v)
{
	int64_t i, m_a = v->n * 2, n_a = 0;
	int32_t max_len = 0;
	mb_anchor_t *a;
	uint8_t *tseq, *qseq2[2];

	for (i = 0; i < v->n; ++i)
		if (max_len < v->a[i].len)
			max_len = v->a[i].len;
	tseq = Kmalloc(km, uint8_t, max_len + qlen * 2);
	qseq2[0] = tseq + max_len;
	qseq2[1] = qseq2[0] + qlen;
	memcpy(qseq2[0], qseq0, qlen);
	for (i = 0; i < qlen; ++i)
		qseq2[1][qlen - 1 - i] = qseq0[i] > 3? 4 : 3 - qseq0[i];

	a = Kmalloc(km, mb_anchor_t, m_a);
	for (i = 0; i < v->n; ++i) {
		mb_anchor_t *p, *q = &v->a[i];
		const l2b_ctg_t *ctg = &l2b->ctg[q->sid>>1];
		int32_t rev = q->sid&1;
		int64_t tpos = q->tpos - (ctg->off * 2 + ctg->len * rev); // NB: requiring concatenated ::tpos!!
		int64_t ts = tpos + 1 - q->len;
		int32_t qs = q->qpos + 1 - q->len, j, j0;
		uint8_t t_allow, q_allow;
		l2b_meth_t mt;
		const uint8_t *qseq = qseq2[rev] + qs;
		l2b_getseq(l2b, q->sid>>1, ts, ts + q->len, tseq);
		mt = q->sid&1? l2b_meth_rev(mt0) : mt0;
		t_allow = mt == L2B_METH_C2T? 1 : 2;
		q_allow = mt == L2B_METH_C2T? 3 : 0;
		for (j0 = j = 0; j <= q->len; ++j) {
			if (j == q->len || tseq[j] == 4 || qseq[j] == 4 || (tseq[j] != qseq[j] && !(tseq[j] == t_allow && qseq[j] == q_allow))) {
				if (j - j0 >= min_len) {
					Kgrow(km, mb_anchor_t, a, n_a, m_a);
					p = &a[n_a++];
					*p = *q;
					p->len = j - j0;
					p->qpos = q->qpos - (q->len - j);
					p->tpos = q->tpos - (q->len - j);
				}
				j0 = j + 1;
			}
		}
	}
	kfree(km, tseq);
	Kgrow(km, mb_anchor_t, v->a, n_a, v->m);
	memcpy(v->a, a, n_a * sizeof(mb_anchor_t));
	v->n = n_a;
	kfree(km, a);
}

/* Converting seed intervals to anchors, given u sorted and deduplicated by
 * anchor_rows() and *sa the SA values of the rows it collected for u. Moves
 * *sa past them, to the next read's. */
double mb_anchor_sa(void *km, const mb_idx_t *idx, const mb_sai_v *u, uint64_t **sa, int32_t min_len, int32_t qlen, const uint8_t *qseq, l2b_meth_t mt, int32_t max_occ, mb_anchor_v *v)
{
	int64_t i, i0, k, n;
	double seed_ratio = 1.0;

	v->n = 0;
	if (u->n == 0) return seed_ratio; // no anchors

	for (i = 0, k = 0; i < u->n; ++i) // pre-calculate the size of v->a
		k += u->a[i].size < max_occ? u->a[i].size : max_occ;
	Kgrow(km, mb_anchor_t, v->a, k - 1, v->m); // preallocate

	for (i = 1, i0 = 0; i <= u->n; ++i) { // intervals sharing ::{x[0],size} share their SA values
		if (i == u->n || u->a[i].x[0] != u->a[i0].x[0] || u->a[i].size != u->a[i0].size) {
			n = anchor_rows1(&u->a[i0], max_occ, 0);
			for (k = 0; k < n; ++k)
				anchor_add(km, idx, u, i0, i, *(*sa)++, qlen, mt, v);
			i0 = i;
		}
	}

	if (mt != L2B_METH_NONE && v->n > 0) {
		int64_t t0, t1;
		for (i = 0, t0 = 0; i < v->n; ++i) t0 += v->a[i].len;
		mb_anchor_split_meth(km, idx->l2b, min_len, qlen, qseq, mt, v);
		for (i = 0, t1 = 0; i < v->n; ++i) t1 += v->a[i].len;
		seed_ratio = (double)t1 / t0;
	}

	radix_sort_mb_anchor(v->a, v->a + v->n);
	for (i = 0; i < v->n; ++i) { // adjust mb_anchor_t::tpos
		mb_anchor_t *q = &v->a[i];
		const l2b_ctg_t *ctg = &idx->l2b->ctg[q->sid>>1];
		q->tpos -= ctg->off * 2 + ctg->len * (q->sid&1);
	}
	mb_anchor_dedup(v);
	return seed_ratio;
}

/* Converting seed intervals to anchors. This function resolves all SA rows
 * of the read together in mb_bwt_sa_batch(), which prefetches the BWT blocks
 * of every pending row before it walks any of them. */
double mb_anchor(void *km, const mb_idx_t *idx, mb_sai_v *u, int32_t min_len, int32_t qlen, const uint8_t *qseq, l2b_meth_t mt, int32_t max_occ, mb_anchor_v *v)
{
	int64_t n_sa = 0, m_sa = 0;
	uint64_t *sa = 0, *p;
	double seed_ratio;
	v->n = 0;
	if (u->n == 0) return 1.0; // no anchors
	anchor_rows(km, u, max_occ, &n_sa, &m_sa, &sa);
	anchor_resolve(km, idx->bwt, n_sa, sa);
	p = sa;
	seed_ratio = mb_anchor_sa(km, idx, u, &p, min_len, qlen, qseq, mt, max_occ, v);
	kfree(km, sa);
	return seed_ratio;
}

/* mb_seed_intv_batch(), then the SA rows of all n_seq reads' anchors resolved
 * together, so that mb_bwt_sa_batch() has many rows in flight. Sorts and
 * deduplicates each v[i] as mb_anchor() does and returns the SA values, read
 * after read, for mb_anchor_sa(); returns NULL when there are none or under
 * --dbg-seed, which prints each read's intervals in their seeding order. */
uint64_t *mb_seed_sa_batch(void *km, const mb_idx_t *idx, int32_t n_seq, const int32_t *len, uint8_t *const* seq, int32_t min_len, int32_t max_sub_occ, int32_t max_occ, mb_sai_v *v)
{
	const mb_bwt_t *bwt = idx->bwt;
	const l2b_t *l2b = idx->l2b;
	int64_t n_sa = 0, m_sa = 0;
	uint64_t *sa = 0;
	const uint64_t *pac = 0;
	int32_t i;
	if (!idx->is_meth && l2b && l2b->pac && bwt->sa && bwt->seq_len == l2b->tot_len * 2 && l2b->n_pac * 32 >= l2b->tot_len)
		pac = l2b->pac; // the BWT text is the l2b's bases followed by their reverse complement
	seed_intv_batch(km, bwt, pac, n_seq, len, seq, min_len, max_sub_occ, v);
	if (kom_dbg_flag & MB_DBG_SEED) return 0;
	for (i = 0; i < n_seq; ++i)
		anchor_rows(km, &v[i], max_occ, &n_sa, &m_sa, &sa);
	anchor_resolve(km, bwt, n_sa, sa);
	return sa;
}

void mb_anchor_sort(const l2b_t *l2b, int64_t n_a, mb_anchor_t *a)
{
	int64_t i;
	if (n_a <= 1) return;
	for (i = 0; i < n_a; ++i) {
		const l2b_ctg_t *ctg = &l2b->ctg[a[i].sid>>1];
		a[i].tpos += ctg->off * 2 + ctg->len * (a[i].sid&1);
	}
	radix_sort_mb_anchor(a, a + n_a);
	for (i = 0; i < n_a; ++i) {
		const l2b_ctg_t *ctg = &l2b->ctg[a[i].sid>>1];
		a[i].tpos -= ctg->off * 2 + ctg->len * (a[i].sid&1);
	}
}
