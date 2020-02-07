#ifndef __BML_THRESHOLD_ELLSORT_H
#define __BML_THRESHOLD_ELLSORT_H

#include "bml_types_ellsort.h"

bml_matrix_ellsort_t *bml_threshold_new_ellsort(
    bml_matrix_ellsort_t * A,
    double threshold);

bml_matrix_ellsort_t
    * bml_threshold_new_ellsort_single_real(bml_matrix_ellsort_t * A,
                                            double threshold);

bml_matrix_ellsort_t
    * bml_threshold_new_ellsort_double_real(bml_matrix_ellsort_t * A,
                                            double threshold);

bml_matrix_ellsort_t
    * bml_threshold_new_ellsort_single_complex(bml_matrix_ellsort_t * A,
                                               double threshold);

bml_matrix_ellsort_t
    * bml_threshold_new_ellsort_double_complex(bml_matrix_ellsort_t * A,
                                               double threshold);

void bml_threshold_ellsort(
    bml_matrix_ellsort_t * A,
    double threshold);

void bml_threshold_ellsort_single_real(
    bml_matrix_ellsort_t * A,
    double threshold);

void bml_threshold_ellsort_double_real(
    bml_matrix_ellsort_t * A,
    double threshold);

void bml_threshold_ellsort_single_complex(
    bml_matrix_ellsort_t * A,
    double threshold);

void bml_threshold_ellsort_double_complex(
    bml_matrix_ellsort_t * A,
    double threshold);

#endif
