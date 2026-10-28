#ifndef MB_BSEQ_H
#define MB_BSEQ_H

#include <stdint.h>
#include <string.h>
#include "kommon.h"

#ifdef __cplusplus
extern "C" {
#endif

struct mb_bseq_file_s;
typedef struct mb_bseq_file_s mb_bseq_file_t;

typedef struct {
	uint64_t l_seq, id; // FIXME: bseq doesn't fully support 64-bit integers yet
	char *name, *seq, *qual, *comment;
} mb_bseq1_t;

mb_bseq_file_t *mb_bseq_open(const char *fn);
void mb_bseq_close(mb_bseq_file_t *fp);
mb_bseq1_t *mb_bseq_read(mb_bseq_file_t *fp, int64_t chunk_size, int with_qual, int with_comment, int frag_mode, int min_cnt, int64_t max_chunk_size, int *n_);
mb_bseq1_t *mb_bseq_read_frag(int n_fp, mb_bseq_file_t **fp, int64_t chunk_size, int with_qual, int with_comment, int *n_);
int mb_bseq_eof(mb_bseq_file_t *fp);
void mb_bseq_warn_buffer(kstring_t *buf); // 0: print warnings to stderr; else append them to buf
void mb_bseq_ahead(mb_bseq_file_t *fp, int on); // 1: inflate ahead of the parser on a thread (regular files only); 0: stop it
void mb_bseq_read_hook(void (*hook)(void *data, const mb_bseq1_t *a, int n), void *data); // 0: none; else hook(data, a, n) after each record a batch takes, a[0..n-1] being the batch so far

static inline int mb_qname_len(const char *s)
{
	int l;
	l = strlen(s);
	return l >= 3 && s[l-1] >= '0' && s[l-1] <= '9' && s[l-2] == '/'? l - 2 : l;
}

static inline int mb_qname_same(const char *s1, const char *s2)
{
	int l1, l2;
	l1 = mb_qname_len(s1);
	l2 = mb_qname_len(s2);
	return (l1 == l2 && strncmp(s1, s2, l1) == 0);
}

#ifdef __cplusplus
}
#endif

#endif
