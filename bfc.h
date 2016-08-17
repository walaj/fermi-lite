#ifndef AC_BFC_H__
#define AC_BFC_H__

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include "htab.h"
#include "kmer.h"
#include "internal.h"
#include "fml.h"
#include "khash.h"

#define _cnt_eq(a, b) ((a)>>14 == (b)>>14)
#define _cnt_hash(a) ((a)>>14)
KHASH_INIT(cnt, uint64_t, char, 0, _cnt_hash, _cnt_eq)
typedef khash_t(cnt) cnthash_t;

struct bfc_ch_s {
  int k;
  cnthash_t **h;
  // private
  int l_pre;
};

typedef struct {
	int n_threads, q, k, l_pre;
	int min_cov; // a k-mer is considered solid if the count is no less than this

	int max_end_ext;
	int win_multi_ec;
	float min_trim_frac;

	// these ec options cannot be changed on the command line
	int w_ec, w_ec_high, w_absent, w_absent_high;
	int max_path_diff, max_heap;
} bfc_opt_t;

/**********************
 *** K-mer counting ***
 **********************/

#define CNT_BUF_SIZE 256

typedef struct { // cache to reduce locking
	uint64_t y[2];
	int is_high;
} insbuf_t;

typedef struct {
	int k, q;
	int n_seqs;
	const fseq1_t *seqs;
	bfc_ch_t *ch;
	int *n_buf;
	insbuf_t **buf;
} cnt_step_t;

float fml_correct_core(const fml_opt_t *opt, int flt_uniq, int n, fseq1_t *seq);

#endif
