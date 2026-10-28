#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "mbpriv.h"
#include "kommon.h"

static char mb_rg_id[256];

/*******************
 * Direct appends  *
 *******************/

/* The per-record writers below append each field directly instead of calling
 * kom_sprintf_lite(), which parses its format string at run time. Every helper
 * keeps room for a trailing NUL, as kom_sprintf_lite() does; the integer helpers
 * print what kom_sprintf_lite()'s %d, %u and %ld print. */

static __attribute__((noinline)) void str_grow(kstring_t *s, size_t l) // str_enlarge()'s rare path, out of line
{
	s->m = s->l + l + 1;
	kom_roundup64(s->m);
	s->s = kom_realloc(char, s->s, s->m);
}

static inline void str_enlarge(kstring_t *s, size_t l)
{
	if (s->l + l + 1 > s->m) str_grow(s, l);
}

static inline void str_putsn(kstring_t *s, const char *p, size_t l)
{
	str_enlarge(s, l);
	memcpy(&s->s[s->l], p, l);
	s->l += l;
}

#define str_putlit(s, lit) str_putsn((s), (lit), sizeof(lit) - 1)

static inline void str_puts(kstring_t *s, const char *p)
{
	str_putsn(s, p, strlen(p));
}

static inline void str_putc(kstring_t *s, int c)
{
	str_enlarge(s, 1);
	s->s[s->l++] = c;
}

static const char mb_digits2[] =
	"00010203040506070809101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899";

static __attribute__((noinline)) void str_putu32(kstring_t *s, uint32_t x) // %u; out of line, so the many call sites stay small
{
	char buf[12], *p = buf + 12;
	while (x >= 100) {
		uint32_t q = x / 100;
		p -= 2, memcpy(p, &mb_digits2[(x - q * 100) * 2], 2);
		x = q;
	}
	if (x >= 10) p -= 2, memcpy(p, &mb_digits2[x * 2], 2);
	else *--p = '0' + x;
	str_putsn(s, p, buf + 12 - p);
}

static __attribute__((noinline)) void str_putu64(kstring_t *s, uint64_t x) // %lu; out of line like str_putu32()
{
	char buf[24], *p = buf + 24;
	while (x >= 100) {
		uint64_t q = x / 100;
		p -= 2, memcpy(p, &mb_digits2[(x - q * 100) * 2], 2);
		x = q;
	}
	if (x >= 10) p -= 2, memcpy(p, &mb_digits2[x * 2], 2);
	else *--p = '0' + x;
	str_putsn(s, p, buf + 24 - p);
}

static inline void str_putd(kstring_t *s, int32_t x) // %d
{
	if (x < 0) str_putc(s, '-'), str_putu32(s, 0U - (uint32_t)x);
	else str_putu32(s, (uint32_t)x);
}

static inline void str_putld(kstring_t *s, long x) // %ld
{
	if (x < 0) str_putc(s, '-'), str_putu64(s, 0UL - (unsigned long)x);
	else str_putu64(s, (unsigned long)x);
}

static inline void str_put_tabd(kstring_t *s, int32_t x) // "\t%d"
{
	str_putc(s, '\t');
	str_putd(s, x);
}

/**************
 * PAF output *
 **************/

static inline void write_tags(kstring_t *s, const mb_hit_t *p)
{
	int32_t nm = p->blen - p->mlen + p->p->n_ambi;
	str_putlit(s, "\tNM:i:"), str_putd(s, nm);
	str_putlit(s, "\tAS:i:"), str_putd(s, p->p->dp_score);
	str_putlit(s, "\tms:i:"), str_putd(s, p->p->dp_max0);
	str_putlit(s, "\tmd:i:"), str_putd(s, p->p->dp_max - p->p->dp_max2);
}

void mb_fmt_paf(kstring_t *s, const l2b_t *l2b, const mb_bseq1_t *t, const mb_hit_t *p, uint64_t opt_flag, int n_seg, int seg_idx)
{
	str_puts(s, t->name);
	if (n_seg > 1 && seg_idx >= 0)
		str_putc(s, '/'), str_putd(s, seg_idx + 1);
	str_putc(s, '\t'), str_putld(s, (long)t->l_seq);
	if (p == 0) { // for unmapped reads
		str_putlit(s, "\t*\t*\t*\t*\t*\t*\t*\t0\t0\t0\n");
		s->s[s->l] = 0;
		return;
	}
	str_put_tabd(s, p->qs);
	str_put_tabd(s, p->qe);
	str_putc(s, '\t'), str_putc(s, p->rev? '-' : '+');
	str_putc(s, '\t'), str_puts(s, l2b->ctg[p->tid].name);
	str_putc(s, '\t'), str_putld(s, (long)l2b->ctg[p->tid].len);
	str_putc(s, '\t'), str_putld(s, (long)p->ts);
	str_putc(s, '\t'), str_putld(s, (long)p->te);
	str_put_tabd(s, p->mlen);
	str_put_tabd(s, p->blen);
	str_put_tabd(s, p->mapq);
	str_putlit(s, "\ttp:A:"), str_putc(s, p->parent == p->id? 'P' : 'S');
	str_putlit(s, "\ts1:i:"), str_putd(s, p->score);
	str_putlit(s, "\tcm:i:"), str_putd(s, p->cnt);
	if (p->parent == p->id) str_putlit(s, "\ts2:i:"), str_putd(s, p->subsc >= 0? p->subsc : 0);
	if (p->p) {
		write_tags(s, p);
		if (p->p->n_cigar > 0) {
			int32_t i;
			str_putlit(s, "\tcg:Z:");
			for (i = 0; i < p->p->n_cigar; ++i)
				str_putu32(s, p->p->cigar[i]>>4), str_putc(s, MB_CIGAR_STR[p->p->cigar[i]&0xf]);
		}
		if (p->p->cs) str_putc(s, '\t'), str_puts(s, (char*)&p->p->cigar[p->p->n_cigar]);
	}
	if ((opt_flag & MB_F_COPY_COMMENT) && t->comment)
		str_putc(s, '\t'), str_puts(s, t->comment);
	str_putc(s, '\n');
	s->s[s->l] = 0;
}

/**************
 * SAM header *
 **************/

char *mb_escape(char *s)
{
	char *p, *q;
	for (p = q = s; *p; ++p) {
		if (*p == '\\') {
			++p;
			if (*p == 't') *q++ = '\t';
			else if (*p == 'n') *q++ = '\n';
			else if (*p == 'r') *q++ = '\r';
			else if (*p == '\\') *q++ = '\\';
			else if (*p == '\0') break;
		} else *q++ = *p;
	}
	*q = '\0';
	return s;
}

static int sam_write_rg_line(kstring_t *str, const char *s)
{
	char *p, *q, *r, *rg_line = 0;
	memset(mb_rg_id, 0, 256);
	if (s == 0) return 0;
	if (strstr(s, "@RG") != s) {
		if (kom_verbose >= 1) fprintf(stderr, "[ERROR] the read group line is not started with @RG\n");
		goto err_set_rg;
	}
	if (strstr(s, "\t") != NULL) {
		if (kom_verbose >= 1) fprintf(stderr, "[ERROR] the read group line contained literal <tab> characters -- replace with escaped tabs: \\t\n");
		goto err_set_rg;
	}
	rg_line = kom_strdup(s);
	mb_escape(rg_line);
	if ((p = strstr(rg_line, "\tID:")) == 0) {
		if (kom_verbose >= 1) fprintf(stderr, "[ERROR] no ID within the read group line\n");
		goto err_set_rg;
	}
	p += 4;
	for (q = p; *q && *q != '\t' && *q != '\n'; ++q);
	if (q - p + 1 > 256) {
		if (kom_verbose >= 1) fprintf(stderr, "[ERROR] @RG:ID is longer than 255 characters\n");
		goto err_set_rg;
	}
	for (q = p, r = mb_rg_id; *q && *q != '\t' && *q != '\n'; ++q)
		*r++ = *q;
	kom_sprintf_lite(str, "%s\n", rg_line);
	free(rg_line);
	return 0;

err_set_rg:
	free(rg_line);
	return -1;
}

int mb_fmt_sam_hdr(kstring_t *str, const l2b_t *idx, const char *rg, const char *ver, int argc, char *argv[])
{
	int i, ret = 0;
	str->l = 0;
	kom_sprintf_lite(str, "@HD\tVN:1.6\tSO:unsorted\tGO:query\n");
	if (idx)
		for (i = 0; i < idx->n_ctg; ++i)
			kom_sprintf_lite(str, "@SQ\tSN:%s\tLN:%ld\n", idx->ctg[i].name, idx->ctg[i].len);
	if (rg) ret = sam_write_rg_line(str, rg);
	kom_sprintf_lite(str, "@PG\tID:minibwa\tPN:minibwa");
	if (ver) kom_sprintf_lite(str, "\tVN:%s", ver);
	if (argc > 1) {
		kom_sprintf_lite(str, "\tCL:minibwa");
		for (i = 0; i < argc; ++i)
			kom_sprintf_lite(str, " %s", argv[i]);
	}
	kom_sprintf_lite(str, "\n");
	return ret;
}

/**************
 * SAM output *
 **************/

static inline void str_copy(kstring_t *s, const char *st, const char *en)
{
	str_enlarge(s, en - st);
	memcpy(&s->s[s->l], st, en - st);
	s->l += en - st;
}

static void sam_write_sq(kstring_t *s, char *seq, int l, int rev, int comp)
{
	if (rev) {
		int i;
		str_enlarge(s, l);
		for (i = 0; i < l; ++i) {
			int c = seq[l - 1 - i];
			s->s[s->l + i] = c < 128 && comp? kom_comp_table[c] : c;
		}
		s->l += l;
	} else str_copy(s, seq, seq + l);
}

static inline const mb_hit_t *get_sam_pri(int n_hit, const mb_hit_t *hit)
{
	int i;
	for (i = 0; i < n_hit; ++i)
		if (hit[i].sam_pri)
			return &hit[i];
	assert(n_hit == 0);
	return NULL;
}

static void write_sam_cigar(kstring_t *s, int sam_flag, int in_tag, int qlen, const mb_hit_t *r, int64_t opt_flag)
{
	if (r->p == 0) {
		str_putc(s, '*');
	} else {
		uint32_t k, clip_len[2];
		clip_len[0] = r->rev? qlen - r->qe : r->qs;
		clip_len[1] = r->rev? r->qs : qlen - r->qe;
		if (in_tag) {
			int clip_char = (((sam_flag&0x800) || ((sam_flag&0x100) && (opt_flag&MB_F_2ND_SEQ))) &&
							 !(opt_flag&MB_F_SUPP_SOFT)) ? 5 : 4;
			str_putlit(s, "\tCG:B:I");
			if (clip_len[0]) str_putc(s, ','), str_putu32(s, clip_len[0]<<4|clip_char);
			for (k = 0; k < r->p->n_cigar; ++k)
				str_putc(s, ','), str_putu32(s, r->p->cigar[k]);
			if (clip_len[1]) str_putc(s, ','), str_putu32(s, clip_len[1]<<4|clip_char);
		} else {
			int clip_char = (((sam_flag&0x800) || ((sam_flag&0x100) && (opt_flag&MB_F_2ND_SEQ))) &&
							 !(opt_flag&MB_F_SUPP_SOFT)) ? 'H' : 'S';
			assert(clip_len[0] < qlen && clip_len[1] < qlen);
			if (clip_len[0]) str_putd(s, (int32_t)clip_len[0]), str_putc(s, clip_char);
			for (k = 0; k < r->p->n_cigar; ++k)
				str_putu32(s, r->p->cigar[k]>>4), str_putc(s, MB_CIGAR_STR[r->p->cigar[k]&0xf]);
			if (clip_len[1]) str_putd(s, (int32_t)clip_len[1]), str_putc(s, clip_char);
		}
	}
}

void mb_fmt_sam(void *km, kstring_t *s, const l2b_t *l2b, const mb_bseq1_t *t, int32_t n_seg, const int32_t *n_hit, mb_hit_t *const*hit, int32_t hit_idx, const mb_opt_t *opt, int seg_idx, int32_t mate_qlen)
{
	int flag, n_h = n_hit[seg_idx];
	int this_tid = -1, this_pos = -1;
	const mb_hit_t *h = hit[seg_idx], *r_prev = NULL, *r_next;
	const mb_hit_t *r = n_h > 0 && hit_idx < n_h && hit_idx >= 0? &h[hit_idx] : NULL;

	assert(n_seg == 1 || n_seg == 2);

	// find the primary of the previous and the next segments, if they are mapped
	if (n_seg > 1) {
		int next_sid = (seg_idx + 1) % n_seg;
		r_prev = r_next = get_sam_pri(n_hit[next_sid], hit[next_sid]);
	} else r_prev = r_next = NULL;

	// write QNAME and FLAG
	str_puts(s, t->name);
	flag = n_seg > 1? 0x1 : 0x0;
	if (r == 0) {
		flag |= 0x4;
	} else {
		if (r->rev) flag |= 0x10;
		if (r->parent != r->id) flag |= 0x100;
		else if (!r->sam_pri) flag |= 0x800;
	}
	if (n_seg > 1) {
		if (r && r->proper_pair) flag |= 0x2;
		if (seg_idx == 0) flag |= 0x40;
		else if (seg_idx == n_seg - 1) flag |= 0x80;
		if (r_next == NULL) flag |= 0x8;
		else if (r_next->rev) flag |= 0x20;
	}
	str_put_tabd(s, flag);

	// write coordinate, MAPQ and CIGAR
	if (r == 0) {
		if (r_prev) {
			this_tid = r_prev->tid, this_pos = r_prev->ts;
			str_putc(s, '\t'), str_puts(s, l2b->ctg[this_tid].name);
			str_put_tabd(s, this_pos+1), str_putlit(s, "\t0\t*");
		} else str_putlit(s, "\t*\t0\t0\t*");
	} else {
		this_tid = r->tid, this_pos = r->ts;
		str_putc(s, '\t'), str_puts(s, l2b->ctg[r->tid].name);
		str_put_tabd(s, (int32_t)(r->ts+1)), str_put_tabd(s, r->mapq), str_putc(s, '\t');
		write_sam_cigar(s, flag, 0, t->l_seq, r, opt->flag);
	}

	// write mate positions
	if (n_seg > 1) {
		int tlen = 0;
		if (this_tid >= 0 && r_next) {
			if (this_tid == r_next->tid) {
				if (r) {
					int this_pos5 = r->rev? r->te - 1 : this_pos;
					int next_pos5 = r_next->rev? r_next->te - 1 : r_next->ts;
					tlen = next_pos5 - this_pos5;
				}
				str_putlit(s, "\t=\t");
			} else str_putc(s, '\t'), str_puts(s, l2b->ctg[r_next->tid].name), str_putc(s, '\t');
			str_putd(s, (int32_t)(r_next->ts + 1)), str_putc(s, '\t');
		} else if (r_next) { // && this_tid < 0
			str_putc(s, '\t'), str_puts(s, l2b->ctg[r_next->tid].name);
			str_put_tabd(s, (int32_t)(r_next->ts + 1)), str_putc(s, '\t');
		} else if (this_tid >= 0) { // && r_next == NULL. In this case, mate is unmapped and placed at this segment's primary alignment
			const mb_hit_t *r_pri;
			r_pri = get_sam_pri(n_h, h);
			if (r->tid == r_pri->tid)
				str_putlit(s, "\t=\t"), str_putd(s, (int32_t)(r_pri->ts + 1)), str_putc(s, '\t');
			else
				str_putc(s, '\t'), str_puts(s, l2b->ctg[r_pri->tid].name), str_put_tabd(s, (int32_t)(r_pri->ts + 1)), str_putc(s, '\t');
		} else str_putlit(s, "\t*\t0\t"); // neither has coordinates
		if (tlen > 0) ++tlen;
		else if (tlen < 0) --tlen;
		str_putd(s, tlen), str_putc(s, '\t');
	} else str_putlit(s, "\t*\t0\t0\t");

	// write SEQ and QUAL
	if (r == 0) {
		sam_write_sq(s, t->seq, t->l_seq, 0, 0);
		str_putc(s, '\t');
		if (t->qual) sam_write_sq(s, t->qual, t->l_seq, 0, 0);
		else str_putc(s, '*');
	} else {
		if ((flag & 0x900) == 0 || (opt->flag & MB_F_SUPP_SOFT)) {
			sam_write_sq(s, t->seq, t->l_seq, r->rev, r->rev);
			str_putc(s, '\t');
			if (t->qual) sam_write_sq(s, t->qual, t->l_seq, r->rev, 0);
			else str_putc(s, '*');
		} else if ((flag & 0x100) && !(opt->flag & MB_F_2ND_SEQ)){
			str_putlit(s, "*\t*");
		} else {
			sam_write_sq(s, t->seq + r->qs, r->qe - r->qs, r->rev, r->rev);
			str_putc(s, '\t');
			if (t->qual) sam_write_sq(s, t->qual + r->qs, r->qe - r->qs, r->rev, 0);
			else str_putc(s, '*');
		}
	}

	// write tags
	if (mb_rg_id[0]) str_putlit(s, "\tRG:Z:"), str_puts(s, mb_rg_id);
	if (n_seg > 2) str_putlit(s, "\tFI:i:"), str_putd(s, seg_idx);
	if (r) {
		write_tags(s, r);
		// MC:Z mate CIGAR and MQ:i mate MAPQ; r_next is the mate's primary (see above).
		if (n_seg > 1 && r_next && r_next->p && r_next->p->n_cigar > 0 && mate_qlen > 0) {
			str_putlit(s, "\tMC:Z:");
			write_sam_cigar(s, 0, 0, mate_qlen, r_next, opt->flag);
			str_putlit(s, "\tMQ:i:"), str_putd(s, r_next->mapq);
		}
		if (r->p->cs) str_putc(s, '\t'), str_puts(s, (char*)&r->p->cigar[r->p->n_cigar]);
		if (r->parent == r->id && r->p && n_h > 1 && h && r >= h && r - h < n_h) { // supplementary aln may exist
			int i, n_sa = 0; // n_sa: number of SA fields
			for (i = 0; i < n_h; ++i)
				if (i != r - h && h[i].parent == h[i].id && h[i].p)
					++n_sa;
			if (n_sa > 0) {
				str_putlit(s, "\tSA:Z:");
				for (i = 0; i < n_h; ++i) {
					const mb_hit_t *q = &h[i];
					int l_M, l_I = 0, l_D = 0, clip5 = 0, clip3 = 0;
					if (r == q || q->parent != q->id || q->p == 0) continue;
					if (q->qe - q->qs < q->te - q->ts) l_M = q->qe - q->qs, l_D = (q->te - q->ts) - l_M;
					else l_M = q->te - q->ts, l_I = (q->qe - q->qs) - l_M;
					clip5 = q->rev? t->l_seq - q->qe : q->qs;
					clip3 = q->rev? q->qs : t->l_seq - q->qe;
					str_puts(s, l2b->ctg[q->tid].name), str_putc(s, ',');
					str_putd(s, (int32_t)(q->ts+1)), str_putc(s, ','), str_putc(s, "+-"[q->rev]), str_putc(s, ',');
					if (clip5) str_putd(s, clip5), str_putc(s, 'S');
					if (l_M) str_putd(s, l_M), str_putc(s, 'M');
					if (l_I) str_putd(s, l_I), str_putc(s, 'I');
					if (l_D) str_putd(s, l_D), str_putc(s, 'D');
					if (clip3) str_putd(s, clip3), str_putc(s, 'S');
					str_putc(s, ','), str_putd(s, q->mapq);
					str_putc(s, ','), str_putd(s, q->blen - q->mlen + q->p->n_ambi), str_putc(s, ';');
				}
			}
			if (opt->xa_max > 0) {
				int i, n_xa = 0;
				for (i = 0; i < n_h; ++i)
					if (i != r - h && h[i].parent == r - h && h[i].p->dp_max >= (double)opt->out_s * r->p->dp_max)
						++n_xa;
				if (n_xa > 0) str_putlit(s, "\tn2:i:"), str_putd(s, n_xa);
				if (n_xa > 0 && n_xa <= opt->xa_max) {
					str_putlit(s, "\tXA:Z:");
					for (i = 0; i < n_h; ++i) {
						const mb_hit_t *q = &h[i];
						if (i != r - h && q->parent == r - h && q->p->dp_max >= (double)opt->out_s * r->p->dp_max) {
							str_puts(s, l2b->ctg[q->tid].name), str_putc(s, ','), str_putc(s, "+-"[q->rev]);
							str_putd(s, (int32_t)(q->ts+1)), str_putc(s, ',');
							write_sam_cigar(s, 0, 0, t->l_seq, q, opt->flag);
							str_putc(s, ','), str_putd(s, q->blen - q->mlen + q->p->n_ambi), str_putc(s, ';');
						}
					}
				}
			}
		}
	}

	if ((opt->flag & MB_F_COPY_COMMENT) && t->comment)
		str_putc(s, '\t'), str_puts(s, t->comment);
	str_putc(s, '\n');
	s->s[s->l] = 0; // we always have room for an extra byte
}

void mb_format(void *km, kstring_t *s, const l2b_t *l2b, const mb_bseq1_t *t, int32_t n_seg, const int32_t *n_hit, mb_hit_t *const*hit, int32_t hit_idx, const mb_opt_t *opt, int seg_idx, int32_t mate_qlen)
{
	if (!(opt->flag & MB_F_PAF))
		mb_fmt_sam(km, s, l2b, t, n_seg, n_hit, hit, hit_idx, opt, seg_idx, mate_qlen);
	else
		mb_fmt_paf(s, l2b, t, hit_idx >= 0? &hit[seg_idx][hit_idx] : 0, opt->flag, n_seg, seg_idx);
}
