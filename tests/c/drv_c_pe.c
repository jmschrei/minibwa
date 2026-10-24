/* Area C driver: the paired-end branch of mb_map_batch().
 *
 * With MB_F_PE, mb_map_batch() pairs reads 2i and 2i+1 using the fixed
 * insert-size statistics in mb_opt_t (pe_avg, pe_std, pe_hi, pe_lo); a
 * trailing odd read stays single. The case writes one line per read with
 * its SAM-primary hit for pytest to compare with `minibwa map -I`, and checks
 * the documented fields of mb_hit_t: one SAM primary per mapped read,
 * coordinates inside the contig, proper_pair agreeing between mates, and no
 * proper_pair or rescued hit when pairing is off. */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "minibwa.h"
#include "minibwa_test.h"

/* args: <idx_prefix> <reads.fa> <avg> <std> <hi> <lo> <pe|nopair|se> <out.tsv>
 * out.tsv: name, index, n_hit, then for the SAM-primary hit: contig, ts, te,
 * strand, mapq, proper_pair, rescued, CIGAR (no clips), qs, qe; then the
 * number of hits with rescued set and the number with proper_pair set. */
static int case_pe_batch(int argc, char *argv[])
{
	mb_opt_t opt;
	mb_idx_t *idx;
	mbt_seqs_t s;
	mb_hit_t **hits;
	int32_t i, j, *n_hit, mode;
	FILE *fp;
	if (argc < 8) return MBT_FAIL("usage: pe_batch <idx> <reads.fa> <avg> <std> <hi> <lo> <pe|nopair|se> <out.tsv>"), 1;
	mb_opt_init(&opt);
	opt.pe_avg = atoi(argv[2]), opt.pe_std = atoi(argv[3]);
	opt.pe_hi = atoi(argv[4]), opt.pe_lo = atoi(argv[5]);
	mode = strcmp(argv[6], "pe") == 0? 0 : strcmp(argv[6], "nopair") == 0? 1 : 2;
	if (mode == 1) opt.flag |= MB_F_NO_PAIRING;
	if (mode == 2) opt.flag &= ~MB_F_PE;
	MBT_REQUIRE(mode != 0 || (opt.flag & MB_F_PE)); // the default preset is paired-end
	idx = mb_idx_load(argv[0], 0);
	MBT_REQUIRE(idx != 0);
	MBT_REQUIRE(mbt_read_seqs(argv[1], &s) > 1);
	n_hit = (int32_t*)calloc(s.n, sizeof(int32_t));
	hits = mb_map_batch(&opt, idx, s.n, s.len, (const char**)s.seq, n_hit, 0, (const char**)s.name);
	MBT_REQUIRE(hits != 0);
	fp = fopen(argv[7], "w");
	MBT_REQUIRE(fp != 0);
	for (i = 0; i < s.n; ++i) {
		const mb_hit_t *pri = 0;
		int32_t n_pri = 0, n_res = 0, n_pp = 0;
		for (j = 0; j < n_hit[i]; ++j) {
			const mb_hit_t *h = &hits[i][j];
			int64_t len = mb_idx_ctg_len(idx, h->tid);
			MBT_CHECKF(len > 0 && 0 <= h->ts && h->ts < h->te && h->te <= len, "read %s: hit %d at %ld-%ld", s.name[i], j, (long)h->ts, (long)h->te);
			MBT_CHECKF(0 <= h->qs && h->qs < h->qe && h->qe <= s.len[i], "read %s: hit %d query %d-%d", s.name[i], j, h->qs, h->qe);
			MBT_CHECKF(h->p != 0 && h->p->n_cigar > 0, "read %s: hit %d has no CIGAR", s.name[i], j);
			if (h->sam_pri) ++n_pri, pri = h;
			n_res += h->rescued, n_pp += h->proper_pair;
		}
		MBT_CHECKF(n_hit[i] == 0 || n_pri == 1, "read %s: %d SAM-primary hits", s.name[i], n_pri);
		if (mode != 0) MBT_CHECKF(n_res == 0 && n_pp == 0, "read %s: pairing is off but %d rescued, %d proper", s.name[i], n_res, n_pp);
		fprintf(fp, "%s\t%d\t%d", s.name[i], i, n_hit[i]);
		if (pri) {
			fprintf(fp, "\t%s\t%ld\t%ld\t%c\t%d\t%d\t%d\t", mb_idx_ctg_name(idx, pri->tid), (long)pri->ts, (long)pri->te, "+-"[pri->rev], pri->mapq, pri->proper_pair, pri->rescued);
			for (j = 0; j < pri->p->n_cigar; ++j)
				fprintf(fp, "%d%c", pri->p->cigar[j]>>4, MB_CIGAR_STR[pri->p->cigar[j]&0xf]);
			fprintf(fp, "\t%d\t%d", pri->qs, pri->qe);
		} else fprintf(fp, "\t*\t-1\t-1\t*\t0\t0\t0\t*\t-1\t-1");
		fprintf(fp, "\t%d\t%d\n", n_res, n_pp);
	}
	fclose(fp);
	for (i = 0; mode == 0 && i + 1 < s.n; i += 2) { // mates agree on proper_pair
		int32_t k, pp[2] = {0, 0};
		for (k = 0; k < 2; ++k)
			for (j = 0; j < n_hit[i+k]; ++j)
				if (hits[i+k][j].sam_pri) pp[k] = hits[i+k][j].proper_pair;
		MBT_CHECKF(pp[0] == pp[1], "pair %s: proper_pair %d and %d", s.name[i], pp[0], pp[1]);
	}
	for (i = 0; i < s.n; ++i) {
		for (j = 0; j < n_hit[i]; ++j) free(hits[i][j].p);
		free(hits[i]);
	}
	free(hits);
	free(n_hit);
	mbt_free_seqs(&s);
	mb_idx_destroy(idx);
	return 0;
}

static const mbt_case_t cases[] = {
	{ "pe_batch", case_pe_batch },
	{ 0, 0 }
};

MBT_MAIN(cases)
