#ifndef MINIBWA_TEST_H
#define MINIBWA_TEST_H

/* A tiny check-and-report framework for C test drivers.
 *
 * A driver defines cases as functions `static int name(int argc, char *argv[])`
 * (argv holds the arguments after the case name), lists them in an
 * mbt_case_t array ending with {0, 0}, and ends with MBT_MAIN(that_array).
 *
 *   driver list              print the case names, one per line
 *   driver run CASE [ARGS]   run one case
 *
 * Each failed check prints a TAP-style "not ok N - ..." line with the file,
 * line and expression. Passing checks are silent unless MBT_VERBOSE is set
 * in the environment. A run ends with "# mbt: CASE checks=N failures=M" and
 * exits 0 only if no check failed and the case returned 0. pytest collects
 * this through minibwa_helpers.run_c_case() / assert_c_case().
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <stdint.h>

typedef int (*mbt_func_t)(int argc, char *argv[]);

typedef struct {
	const char *name;
	mbt_func_t func;
} mbt_case_t;

static int mbt_n_check = 0, mbt_n_fail = 0, mbt_verbose = 0;

static int mbt_report(int ok, const char *file, int line, const char *fmt, ...)
{
	va_list ap;
	++mbt_n_check;
	if (!ok) ++mbt_n_fail;
	if (ok && !mbt_verbose) return ok;
	printf("%s %d - %s:%d: ", ok? "ok" : "not ok", mbt_n_check, file, line);
	va_start(ap, fmt);
	vprintf(fmt, ap);
	va_end(ap);
	putchar('\n');
	fflush(stdout);
	return ok;
}

/* each macro evaluates to 1 if the check passed and 0 otherwise */
#define MBT_CHECK(cond) mbt_report(!!(cond), __FILE__, __LINE__, "%s", #cond)
#define MBT_CHECKF(cond, ...) mbt_report(!!(cond), __FILE__, __LINE__, __VA_ARGS__)
#define MBT_CHECK_INT(a, op, b) do { \
		long long mbt_a_ = (long long)(a), mbt_b_ = (long long)(b); \
		mbt_report(mbt_a_ op mbt_b_, __FILE__, __LINE__, "%s %s %s (%lld vs %lld)", #a, #op, #b, mbt_a_, mbt_b_); \
	} while (0)
#define MBT_CHECK_STR(a, b) do { \
		const char *mbt_a_ = (a), *mbt_b_ = (b); \
		int mbt_ok_ = (mbt_a_ && mbt_b_ && strcmp(mbt_a_, mbt_b_) == 0) || (mbt_a_ == 0 && mbt_b_ == 0); \
		mbt_report(mbt_ok_, __FILE__, __LINE__, "%s == %s (\"%s\" vs \"%s\")", #a, #b, mbt_a_? mbt_a_ : "(null)", mbt_b_? mbt_b_ : "(null)"); \
	} while (0)
#define MBT_FAIL(...) mbt_report(0, __FILE__, __LINE__, __VA_ARGS__)
/* abort the case (after reporting) if a precondition does not hold */
#define MBT_REQUIRE(cond) do { if (!MBT_CHECK(cond)) return 1; } while (0)

static void mbt_note(const char *fmt, ...)
{
	va_list ap;
	fputs("# ", stdout);
	va_start(ap, fmt);
	vprintf(fmt, ap);
	va_end(ap);
	putchar('\n');
}

static int mbt_main(int argc, char *argv[], const mbt_case_t *cases)
{
	const mbt_case_t *c;
	const char *v;
	int ret;
	if (argc >= 2 && strcmp(argv[1], "list") == 0) {
		for (c = cases; c->name; ++c) puts(c->name);
		return 0;
	}
	if (argc < 3 || strcmp(argv[1], "run") != 0) {
		fprintf(stderr, "Usage: %s list | run <case> [args...]\n", argv[0]);
		return 2;
	}
	for (c = cases; c->name; ++c)
		if (strcmp(c->name, argv[2]) == 0) break;
	if (c->name == 0) {
		fprintf(stderr, "ERROR: unknown case '%s'\n", argv[2]);
		return 2;
	}
	v = getenv("MBT_VERBOSE");
	mbt_verbose = v && *v && strcmp(v, "0") != 0;
	ret = c->func(argc - 3, argv + 3);
	if (ret != 0) mbt_note("case returned %d", ret);
	printf("# mbt: %s checks=%d failures=%d\n", c->name, mbt_n_check, mbt_n_fail);
	return mbt_n_fail > 0 || ret != 0? 1 : 0;
}

#define MBT_MAIN(cases) int main(int argc, char *argv[]) { return mbt_main(argc, argv, (cases)); }

/* read the sequences of a FASTA/FASTQ file without kseq; for small test inputs only */
typedef struct {
	int32_t n;
	char **name, **seq;
	int32_t *len;
} mbt_seqs_t;

static int mbt_read_seqs(const char *fn, mbt_seqs_t *s)
{
	FILE *fp;
	char *line = 0;
	size_t cap = 0;
	int32_t m = 0, mode = 0; // mode: 1 after a FASTA header, 2 after a FASTQ header, 3 at FASTQ '+', 4 skip quality
	long l;
	memset(s, 0, sizeof(*s));
	if ((fp = fopen(fn, "r")) == 0) return -1;
	while (1) {
		int c;
		size_t n = 0;
		while ((c = fgetc(fp)) != EOF && c != '\n') {
			if (n + 2 > cap) cap = cap? cap * 2 : 256, line = (char*)realloc(line, cap);
			line[n++] = (char)c;
		}
		if (c == EOF && n == 0) break;
		if (line == 0) cap = 256, line = (char*)malloc(cap);
		line[n] = 0;
		l = (long)n;
		if (mode == 4) { mode = 0; continue; }
		if (mode == 3) { mode = 4; continue; }
		if ((line[0] == '>' || line[0] == '@') && (mode == 0 || mode == 1)) {
			char *p;
			if (s->n == m) {
				m = m? m * 2 : 16;
				s->name = (char**)realloc(s->name, m * sizeof(char*));
				s->seq = (char**)realloc(s->seq, m * sizeof(char*));
				s->len = (int32_t*)realloc(s->len, m * sizeof(int32_t));
			}
			for (p = line + 1; *p && *p != ' ' && *p != '\t'; ++p);
			*p = 0;
			s->name[s->n] = strdup(line + 1);
			s->seq[s->n] = strdup("");
			s->len[s->n++] = 0;
			mode = line[0] == '>'? 1 : 2;
		} else if (mode == 1 || mode == 2) {
			int32_t k = s->n - 1;
			s->seq[k] = (char*)realloc(s->seq[k], s->len[k] + l + 1);
			memcpy(s->seq[k] + s->len[k], line, l + 1);
			s->len[k] += (int32_t)l;
			if (mode == 2) mode = 3;
		}
		if (c == EOF) break;
	}
	free(line);
	fclose(fp);
	return s->n;
}

static void mbt_free_seqs(mbt_seqs_t *s)
{
	int32_t i;
	for (i = 0; i < s->n; ++i) free(s->name[i]), free(s->seq[i]);
	free(s->name); free(s->seq); free(s->len);
	memset(s, 0, sizeof(*s));
}

#endif
