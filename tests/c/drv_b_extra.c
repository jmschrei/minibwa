/* Area B driver: base-alignment bookkeeping in align.c and cs.c —
 * mb_update_extra (CIGAR clean-up, mlen/blen/n_ambi, the max-scoring segment,
 * --eqx, and the cs/ds/MD strings it appends), mb_append_cigar,
 * mb_write_cs_ds, mb_write_MD and mb_update_dp_max — on short hand-written
 * alignments whose expected values are worked out in the comments. */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "mbpriv.h"
#include "ksw2.h"
#include "minibwa_test.h"

// defined in align.c without a prototype in a header
void mb_update_dp_max(int qlen, int n_regs, mb_hit_t *regs, double frac, int a, int b);

static void be_enc(const char *s, uint8_t *x)
{
	for (; *s; ++s) *x++ = kom_nt4_table[(uint8_t)*s];
}

static int be_parse_cigar(const char *s, uint32_t *c)
{
	int n = 0;
	while (*s) {
		char *e;
		long l = strtol(s, &e, 10);
		const char *op = strchr(MB_CIGAR_STR, *e);
		c[n++] = (uint32_t)l << 4 | (uint32_t)(op - MB_CIGAR_STR);
		s = e + 1;
	}
	return n;
}

static void be_cigar_str(const mb_hit_t *r, char *out)
{
	uint32_t k;
	*out = 0;
	for (k = 0; k < (uint32_t)r->p->n_cigar; ++k)
		out += sprintf(out, "%u%c", r->p->cigar[k] >> 4, MB_CIGAR_STR[r->p->cigar[k] & 0xf]);
}

/* a hit with the given CIGAR over [0, qspan) x [0, tspan) */
static void be_hit(mb_hit_t *r, const char *cigar, int qlen, int rev)
{
	uint32_t c[64];
	int n = be_parse_cigar(cigar, c), k, qs = 0, ts = 0;
	memset(r, 0, sizeof(*r));
	mb_append_cigar(r, n, c);
	for (k = 0; k < n; ++k) {
		int op = c[k] & 0xf, len = c[k] >> 4;
		if (op == MB_CIGAR_MATCH || op == MB_CIGAR_EQ_MATCH || op == MB_CIGAR_X_MISMATCH) qs += len, ts += len;
		else if (op == MB_CIGAR_INS) qs += len;
		else ts += len;
	}
	r->rev = rev;
	r->ts = 100, r->te = 100 + ts;
	if (!rev) r->qs = 0, r->qe = qs;
	else r->qs = qlen - qs, r->qe = qlen;
}

static const char *be_tag(const mb_hit_t *r) { return r->p->cs? (const char*)&r->p->cigar[r->p->n_cigar] : ""; }

typedef struct {
	const char *q, *t, *cigar;
	uint64_t flag;
	int log_gap;
	const char *cigar_out, *tag;
	int mlen, blen, n_ambi, dp_max0, dq, dt; // dq/dt: shifts of the start coordinates
} be_case_t;

static void be_run(const be_case_t *c, int rev)
{
	mb_hit_t r;
	uint8_t q[256], t[256];
	int8_t mat[25];
	char cig[512];
	int qlen = (int)strlen(c->q), qs0, ts0, qe0;
	be_enc(c->q, q), be_enc(c->t, t);
	ksw_gen_nt4_mat(mat, 2, 8, 0, 1, 0);
	be_hit(&r, c->cigar, qlen, rev);
	qs0 = r.qs, qe0 = r.qe, ts0 = r.ts;
	mb_update_extra(0, &r, q, t, mat, 12, 2, c->flag, c->log_gap);
	be_cigar_str(&r, cig);
	MBT_CHECKF(strcmp(cig, c->cigar_out) == 0, "%s (rev %d): CIGAR %s, expected %s", c->cigar, rev, cig, c->cigar_out);
	MBT_CHECKF(r.mlen == c->mlen && r.blen == c->blen && (int)r.p->n_ambi == c->n_ambi, "%s: mlen %d blen %d n_ambi %d, expected %d %d %d", c->cigar, r.mlen, r.blen, r.p->n_ambi, c->mlen, c->blen, c->n_ambi);
	MBT_CHECKF(r.p->dp_max0 == c->dp_max0 && r.p->dp_max == c->dp_max0, "%s: dp_max0 %d dp_max %d, expected %d", c->cigar, r.p->dp_max0, r.p->dp_max, c->dp_max0);
	MBT_CHECKF(strcmp(be_tag(&r), c->tag) == 0, "%s: tag \"%s\", expected \"%s\"", c->cigar, be_tag(&r), c->tag);
	MBT_CHECKF(r.ts == ts0 + c->dt, "%s: ts moved by %ld, expected %d", c->cigar, (long)(r.ts - ts0), c->dt);
	if (!rev) MBT_CHECKF(r.qs == qs0 + c->dq && r.qe == qe0, "%s: qs moved by %d, expected %d", c->cigar, r.qs - qs0, c->dq);
	else MBT_CHECKF(r.qe == qe0 - c->dq && r.qs == qs0, "%s (rev): qe moved by %d, expected %d", c->cigar, qe0 - r.qe, c->dq);
	free(r.p);
}

static int case_update_extra(int argc, char *argv[])
{
	static const be_case_t cases[] = {
		// identical: 10 matches, segment score 20
		{ "ACGTACGTAC", "ACGTACGTAC", "10M", 0, 0, "10M", "", 10, 10, 0, 20, 0, 0 },
		// a mismatch at 3 and N at 6: blen 9 (N excluded), mlen 8; segment 6-8 -> 0, then 2+2-1+2+2+2 = 9
		{ "ACGAACNTAC", "ACGTACGTAC", "10M", MB_F_WRITE_CS, 0, "10M", "cs:Z::3*ta:2*gn:3", 8, 9, 1, 9, 0, 0 },
		{ "ACGAACNTAC", "ACGTACGTAC", "10M", MB_F_WRITE_MD, 0, "10M", "MD:Z:3T2G3", 8, 9, 1, 9, 0, 0 },
		// N against N is not a difference in cs or MD, but it is ambiguous (n_ambi 1, not in mlen/blen)
		{ "ACGN", "ACTN", "4M", MB_F_WRITE_CS, 0, "4M", "cs:Z::2*tg:1", 2, 3, 1, 4, 0, 0 },
		{ "ACGN", "ACTN", "4M", MB_F_WRITE_MD, 0, "4M", "MD:Z:2T1", 2, 3, 1, 4, 0, 0 },
		// --eqx: '=' runs and 'X' runs; N counts as X even against N
		{ "ACGN", "ACTN", "4M", MB_F_EQX, 0, "2=2X", "", 2, 3, 1, 4, 0, 0 },
		{ "ACGTTACG", "ACGTTACG", "8M", MB_F_EQX, 0, "8=", "", 8, 8, 0, 16, 0, 0 },
		{ "ACGTCACG", "ACGTACG", "4M1I3M", MB_F_EQX | MB_F_WRITE_CS, 0, "4=1I3=", "cs:Z::4+c:3", 7, 8, 0, 8, 0, 0 },
		{ "ACGACATG", "ACGTCATG", "8M", MB_F_EQX | MB_F_WRITE_MD, 0, "3=1X4=", "MD:Z:3T4", 7, 8, 0, 8, 0, 0 },
		// an insertion in a T run moves to the run's left end: 5M1I4M -> 3M1I6M; ds brackets the repeated base
		{ "ACGTTTTGCA", "ACGTTTGCA", "5M1I4M", MB_F_WRITE_CS, 0, "3M1I6M", "cs:Z::3+t:6", 9, 10, 0, 12, 0, 0 },
		{ "ACGTTTTGCA", "ACGTTTGCA", "5M1I4M", MB_F_WRITE_DS, 0, "3M1I6M", "ds:Z::3+[t]:6", 9, 10, 0, 12, 0, 0 },
		{ "ACGTTTTGCA", "ACGTTTGCA", "5M1I4M", MB_F_WRITE_MD, 1, "3M1I6M", "MD:Z:9", 9, 10, 0, 12, 0, 0 },
		// a 3 bp deletion in CTAGTAGCC: 2M3D4M left-aligns to 1M3D5M (TAG deleted); ds: all of it repeats
		{ "CTAGCC", "CTAGTAGCC", "2M3D4M", MB_F_WRITE_CS, 0, "1M3D5M", "cs:Z::1-tag:5", 6, 9, 0, 10, 0, 0 },
		{ "CTAGCC", "CTAGTAGCC", "2M3D4M", MB_F_WRITE_DS, 0, "1M3D5M", "ds:Z::1-[tag]:5", 6, 9, 0, 10, 0, 0 },
		{ "CTAGCC", "CTAGTAGCC", "2M3D4M", MB_F_WRITE_MD, 0, "1M3D5M", "MD:Z:1^TAG5", 6, 9, 0, 10, 0, 0 },
		// partial repeat context: insertion CGTT after C, followed by CGAA: "[cg]tt"
		{ "CCGTTCGAA", "CCGAA", "2M4I3M", MB_F_WRITE_DS, 0, "1M4I4M", "ds:Z::1+[cg]tt:4", 5, 9, 0, 8, 0, 0 },
		// 5I6D7I-style runs collapse to one I and one D: 3M2I3D1I3M -> 3M3I3D3M
		{ "ACGTTGCAT", "ACGGGGCAT", "3M2I3D1I3M", MB_F_WRITE_CS, 0, "3M3I3D3M", "cs:Z::3+ttg-ggg:3", 6, 12, 0, 6, 0, 0 },
		// a leading insertion or deletion is dropped and the start moves (forward: qs; reverse: qe)
		{ "TTACGTACGT", "ACGTACGT", "2I8M", MB_F_WRITE_CS, 0, "8M", "cs:Z::8", 8, 8, 0, 16, 2, 0 },
		{ "ACGTACGT", "GGACGTACGT", "2D8M", MB_F_WRITE_MD, 0, "8M", "MD:Z:8", 8, 8, 0, 16, 0, 2 },
		// a deletion of CCCC after ...TAC moves left by one (the C before it) and the segment restarts after it:
		// 10 - 14 -> 0, then 7 matches = 14
		{ "ACGTACGTACGT", "ACGTACCCCCGTACGT", "6M4D6M", 0, 0, "5M4D7M", "", 12, 16, 0, 14, 0, 0 },
		// the segment score charges a gap q+e = 14 (short reads) or q + e*log2(1+len) (long reads):
		// 20 - 14 + 20 = 26 vs 20 - 16.6 + 20 = 23.4 -> 23
		{ "ACGTACGTACACGTACGTAC", "ACGTACGTACTTTTACGTACGTAC", "10M4D10M", 0, 0, "10M4D10M", "", 20, 24, 0, 26, 0, 0 },
		{ "ACGTACGTACACGTACGTAC", "ACGTACGTACTTTTACGTACGTAC", "10M4D10M", 0, 1, "10M4D10M", "", 20, 24, 0, 23, 0, 0 },
	};
	int i;
	for (i = 0; i < (int)(sizeof(cases) / sizeof(cases[0])); ++i) {
		be_run(&cases[i], 0);
		if (cases[i].dq || cases[i].dt) be_run(&cases[i], 1);
	}
	return 0;
}

/* mb_append_cigar: allocates on first use, merges equal operations at the junction, grows */
static int case_append_cigar(int argc, char *argv[])
{
	mb_hit_t r;
	uint32_t c[3];
	int i;
	char s[8192];
	memset(&r, 0, sizeof(r));
	mb_append_cigar(&r, 0, c);
	MBT_CHECK(r.p == 0);
	c[0] = 5<<4 | MB_CIGAR_MATCH;
	mb_append_cigar(&r, 1, c);
	MBT_REQUIRE(r.p != 0 && r.p->n_cigar == 1 && r.p->cap >= 1 + sizeof(mb_extra_t)/4);
	c[0] = 3<<4 | MB_CIGAR_MATCH, c[1] = 2<<4 | MB_CIGAR_INS;
	mb_append_cigar(&r, 2, c);
	c[0] = 2<<4 | MB_CIGAR_INS, c[1] = 1<<4 | MB_CIGAR_DEL, c[2] = 4<<4 | MB_CIGAR_MATCH;
	mb_append_cigar(&r, 3, c);
	be_cigar_str(&r, s);
	MBT_CHECK_STR(s, "8M4I1D4M");
	for (i = 0; i < 300; ++i) {
		c[0] = 1<<4 | MB_CIGAR_DEL, c[1] = 1<<4 | MB_CIGAR_MATCH;
		mb_append_cigar(&r, 2, c);
	}
	MBT_CHECK_INT(r.p->n_cigar, ==, 4 + 600);
	MBT_CHECK(r.p->cap >= r.p->n_cigar + sizeof(mb_extra_t)/4);
	MBT_CHECK(r.p->cigar[603] == (1<<4 | MB_CIGAR_MATCH) && r.p->cigar[3] == (4<<4 | MB_CIGAR_MATCH));
	free(r.p);
	return 0;
}

/* mb_write_cs_ds and mb_write_MD called directly, including an intron ('N' op, '~' in cs,
 * skipped by MD) that minibwa itself never produces */
static int case_write_tags(int argc, char *argv[])
{
	mb_hit_t r;
	uint8_t q[64], t[64];
	kstring_t s = {0, 0, 0};
	be_enc("ACGTAGGCA", q);
	be_enc("ACGTGTAAAAAAAAGAGGCA", t);
	be_hit(&r, "4M11N5M", 9, 0);
	mb_write_cs_ds(0, &s, t, q, &r, 0);
	MBT_CHECK_STR(s.s, "cs:Z::4~gt11ag:5");
	s.l = 0;
	mb_write_MD(0, &s, t, q, &r);
	MBT_CHECK_STR(s.s, "MD:Z:9");
	free(r.p);
	// both ends of a deletion in its repeat context: deleting "CA" from ACACAC flanked by G
	be_enc("GACACG", q);
	be_enc("GACACACG", t);
	be_hit(&r, "1M2D5M", 6, 0);
	s.l = 0;
	mb_write_cs_ds(0, &s, t, q, &r, 1);
	MBT_CHECK_STR(s.s, "ds:Z::1-[ac]:5");
	s.l = 0;
	mb_write_cs_ds(0, &s, t, q, &r, 0);
	MBT_CHECK_STR(s.s, "cs:Z::1-ac:5");
	free(r.p);
	// insertion AGCT between CT and AC: its suffix CT repeats the two bases before it, its prefix A
	// the base after it, so "[a]g[ct]"
	be_enc("CTAGCTAC", q);
	be_enc("CTAC", t);
	be_hit(&r, "2M4I2M", 8, 0);
	s.l = 0;
	mb_write_cs_ds(0, &s, t, q, &r, 1);
	MBT_CHECK_STR(s.s, "ds:Z::2+[a]g[ct]:2");
	free(r.p);
	// CAC between GGC and CAT: suffix C repeats the base before, prefix CA the bases after;
	// together they make up the whole insertion, so it is bracketed whole
	be_enc("GGCCACCAT", q);
	be_enc("GGCCAT", t);
	be_hit(&r, "3M3I3M", 9, 0);
	s.l = 0;
	mb_write_cs_ds(0, &s, t, q, &r, 1);
	MBT_CHECK_STR(s.s, "ds:Z::3+[cac]:3");
	free(r.p);
	be_enc("CCGTTCGAA", q);
	be_enc("CCGAA", t);
	be_hit(&r, "2M4I3M", 9, 0); // not left-aligned: suffix c repeats the base before, prefix g the base after
	s.l = 0;
	mb_write_cs_ds(0, &s, t, q, &r, 1);
	MBT_CHECK_STR(s.s, "ds:Z::2+[g]tt[c]:3");
	free(r.p);
	free(s.s);
	return 0;
}

static mb_hit_t be_dp_hit(int qs, int qe, int mlen, int blen, const char *cigar, int dp_max)
{
	mb_hit_t r;
	be_hit(&r, cigar, qe, 0);
	r.qs = qs, r.qe = qe, r.mlen = mlen, r.blen = blen;
	r.p->dp_max = r.p->dp_max0 = dp_max;
	return r;
}

/* mb_update_dp_max (long reads): when the best hit covers >= 90% of the read and the second
 * covers >= sqrt(0.9) of the best, every DP max is recomputed as
 * a * (mlen - b2 * n_mis - sum b2*log2(1+gap)), n_mis counting mismatches plus the
 * clipped bases divided by b2, and b2 = 0.5 / max(0.02, 1 - event identity of the best),
 * replaced by a/b when b2 * a < b. */
static int case_update_dp_max(int argc, char *argv[])
{
	mb_hit_t h[2];
	// identity 1 -> b2 = 25; the second hit has one mismatch: 2 * (399 - 25) = 748
	h[0] = be_dp_hit(0, 400, 400, 400, "400M", 800);
	h[1] = be_dp_hit(0, 400, 399, 400, "400M", 790);
	mb_update_dp_max(400, 2, h, 0.9, 2, 8);
	MBT_CHECK(h[0].p->dp_max == 800 && h[1].p->dp_max == 748);
	MBT_CHECK(h[0].p->dp_max0 == 800 && h[1].p->dp_max0 == 790);
	free(h[0].p), free(h[1].p);
	// a 4 bp insertion: event identity 400/401 -> b2 = 25; 2 * (400 - 25 * log2(5)) = 684
	// second: 14 bases clipped and 10 mismatches: n_mis = 10 + (int)(14/25 + .499) = 11 -> 2 * (380 - 275) = 210
	h[0] = be_dp_hit(0, 404, 400, 404, "199M4I201M", 790);
	h[1] = be_dp_hit(0, 390, 380, 390, "390M", 700);
	mb_update_dp_max(404, 2, h, 0.9, 2, 8);
	MBT_CHECK_INT(h[0].p->dp_max, ==, 684);
	MBT_CHECK_INT(h[1].p->dp_max, ==, 210);
	free(h[0].p), free(h[1].p);
	// low identity (0.75): 0.5/0.25 = 2 gives b2 * a = 4 < b = 8, so b2 = a/b = 0.25
	h[0] = be_dp_hit(0, 400, 300, 400, "400M", 600);
	h[1] = be_dp_hit(0, 400, 280, 400, "400M", 560);
	mb_update_dp_max(400, 2, h, 0.9, 2, 8);
	MBT_CHECK(h[0].p->dp_max == 550 && h[1].p->dp_max == 500);
	free(h[0].p), free(h[1].p);
	// the best hit covers less than 90% of the read: nothing changes
	h[0] = be_dp_hit(0, 350, 350, 350, "350M", 700);
	h[1] = be_dp_hit(0, 350, 349, 350, "350M", 690);
	mb_update_dp_max(400, 2, h, 0.9, 2, 8);
	MBT_CHECK(h[0].p->dp_max == 700 && h[1].p->dp_max == 690);
	// a single hit: nothing changes
	mb_update_dp_max(350, 1, h, 0.9, 2, 8);
	MBT_CHECK(h[0].p->dp_max == 700);
	free(h[0].p), free(h[1].p);
	return 0;
}

static const mbt_case_t cases[] = {
	{ "update_extra", case_update_extra },
	{ "append_cigar", case_append_cigar },
	{ "write_tags", case_write_tags },
	{ "update_dp_max", case_update_dp_max },
	{ 0, 0 }
};

MBT_MAIN(cases)
