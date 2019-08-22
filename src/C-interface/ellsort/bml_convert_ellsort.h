/** \file */

#ifndef __BML_CONVERT_ELLSORT_H
#define __BML_CONVERT_ELLSORT_H

#include "bml_types_ellsort.h"

bml_matrix_ellsort_t *bml_convert_ellsort(
    bml_matrix_t * A,
    bml_matrix_precision_t matrix_precision,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_ellsort_t *bml_convert_ellsort_single_real(
    bml_matrix_t * A,
    bml_matrix_precision_t matrix_precision,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_ellsort_t *bml_convert_ellsort_double_real(
    bml_matrix_t * A,
    bml_matrix_precision_t matrix_precision,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_ellsort_t *bml_convert_ellsort_single_complex(
    bml_matrix_t * A,
    bml_matrix_precision_t matrix_precision,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_ellsort_t *bml_convert_ellsort_double_complex(
    bml_matrix_t * A,
    bml_matrix_precision_t matrix_precision,
    int M,
    bml_distribution_mode_t distrib_mode);

#endif
