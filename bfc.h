#ifndef AC_BFC_H
#define AC_BFC_H

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include "htab.h"
#include "kmer.h"
#include "internal.h"
#include "fml.h"

float fml_correct_core(const fml_opt_t *opt, int flt_uniq, int n, fseq1_t *seq);

#endif
