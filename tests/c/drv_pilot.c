/* Pilot C driver: map simulated reads through the libminibwa API (mb_map) the
 * way api-test/ex-one.c does, and check each read's SAM-primary hit against
 * the simulator's truth. */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "minibwa.h"
#include "minibwa_test.h"

typedef struct {
	char name[256], ctg[256], strand;
	long pos, end;
} truth_t;

static int read_truth(const char *fn, int n_max, truth_t *t)
{
	FILE *fp;
	int n = 0;
	if ((fp = fopen(fn, "r")) == 0) return -1;
	while (n < n_max && fscanf(fp, "%255s %255s %ld %ld %c", t[n].name, t[n].ctg, &t[n].pos, &t[n].end, &t[n].strand) == 5)
		++n;
	fclose(fp);
	return n;
}

/* args: <idx_prefix> <reads.fa> <truth.tsv> <max_dist>
 * truth.tsv: name, contig, 0-based start, end, strand per read */
static int case_map_one(int argc, char *argv[])
{
	mb_opt_t opt;
	mb_idx_t *idx;
	mbt_seqs_t s;
	truth_t *t;
	int32_t i, n_t, max_dist;
	if (argc < 4) return MBT_FAIL("usage: map_one <idx> <reads.fa> <truth.tsv> <max_dist>"), 1;
	max_dist = atoi(argv[3]);
	mb_opt_init(&opt);
	idx = mb_idx_load(argv[0], 0);
	MBT_REQUIRE(idx != 0);
	MBT_REQUIRE(mbt_read_seqs(argv[1], &s) > 0);
	t = (truth_t*)calloc(s.n, sizeof(truth_t));
	n_t = read_truth(argv[2], s.n, t);
	MBT_CHECK_INT(n_t, ==, s.n);
	MBT_CHECK(mb_idx_ctg_name(idx, -1) == 0);
	MBT_CHECK_INT(mb_idx_ctg_len(idx, -1), ==, -1);
	for (i = 0; i < s.n && i < n_t; ++i) {
		mb_hit_t *hit;
		int32_t j, n_hit, n_pri = 0;
		MBT_CHECK_STR(s.name[i], t[i].name);
		hit = mb_map(&opt, idx, s.len[i], s.seq[i], 0, &n_hit, 0, s.name[i]);
		MBT_CHECKF(n_hit > 0, "read %s has no hit", s.name[i]);
		for (j = 0; j < n_hit; ++j) {
			mb_hit_t *h = &hit[j];
			const char *name = mb_idx_ctg_name(idx, h->tid);
			int64_t ctg_len = mb_idx_ctg_len(idx, h->tid);
			long start5;
			MBT_CHECKF(name != 0 && ctg_len > 0, "read %s: hit %d has tid %ld out of range", s.name[i], j, (long)h->tid);
			MBT_CHECKF(0 <= h->ts && h->ts < h->te && h->te <= ctg_len, "read %s: hit %d spans %ld-%ld on a %ld bp contig", s.name[i], j, (long)h->ts, (long)h->te, (long)ctg_len);
			MBT_CHECKF(0 <= h->qs && h->qs < h->qe && h->qe <= s.len[i], "read %s: hit %d query span %d-%d", s.name[i], j, h->qs, h->qe);
			MBT_CHECKF(h->p != 0 && h->p->n_cigar > 0, "read %s: hit %d has no base alignment", s.name[i], j);
			if (!h->sam_pri) continue;
			++n_pri;
			start5 = (long)h->ts - (h->rev? s.len[i] - h->qe : h->qs); // leftmost reference position of the unclipped read
			MBT_CHECK_STR(name, t[i].ctg);
			MBT_CHECKF("+-"[h->rev] == t[i].strand, "read %s: strand %c, truth %c", s.name[i], "+-"[h->rev], t[i].strand);
			MBT_CHECKF(labs(start5 - t[i].pos) <= max_dist, "read %s: unclipped start %ld, truth %ld", s.name[i], start5, t[i].pos);
		}
		MBT_CHECKF(n_hit == 0 || n_pri == 1, "read %s: %d hits flagged primary", s.name[i], n_pri);
		for (j = 0; j < n_hit; ++j) free(hit[j].p);
		free(hit);
	}
	free(t);
	mbt_free_seqs(&s);
	mb_idx_destroy(idx);
	return 0;
}

static const mbt_case_t cases[] = {
	{ "map_one", case_map_one },
	{ 0, 0 }
};

MBT_MAIN(cases)
