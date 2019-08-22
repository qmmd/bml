#ifndef __BML_SCALE_ELLSORT_H
#define __BML_SCALE_ELLSORT_H

#include "bml_types_ellsort.h"

#include <complex.h>

bml_matrix_ellsort_t *bml_scale_ellsort_new(
    void *scale_factor,
    bml_matrix_ellsort_t * A);

bml_matrix_ellsort_t *bml_scale_ellsort_new_single_real(
    void *scale_factor,
    bml_matrix_ellsort_t * A);

bml_matrix_ellsort_t *bml_scale_ellsort_new_double_real(
    void *scale_factor,
    bml_matrix_ellsort_t * A);

bml_matrix_ellsort_t *bml_scale_ellsort_new_single_complex(
    void *scale_factor,
    bml_matrix_ellsort_t * A);

bml_matrix_ellsort_t *bml_scale_ellsort_new_double_complex(
    void *scale_factor,
    bml_matrix_ellsort_t * A);

void bml_scale_ellsort(
    void *scale_factor,
    bml_matrix_ellsort_t * A,
    bml_matrix_ellsort_t * B);

void bml_scale_ellsort_single_real(
    void *scale_factor,
    bml_matrix_ellsort_t * A,
    bml_matrix_ellsort_t * B);

void bml_scale_ellsort_double_real(
    void *scale_factor,
    bml_matrix_ellsort_t * A,
    bml_matrix_ellsort_t * B);

void bml_scale_ellsort_single_complex(
    void *scale_factor,
    bml_matrix_ellsort_t * A,
    bml_matrix_ellsort_t * B);

void bml_scale_ellsort_double_complex(
    void *scale_factor,
    bml_matrix_ellsort_t * A,
    bml_matrix_ellsort_t * B);

void bml_scale_inplace_ellsort(
    void *scale_factor,
    bml_matrix_ellsort_t * A);

void bml_scale_inplace_ellsort_single_real(
    void *scale_factor,
    bml_matrix_ellsort_t * A);

void bml_scale_inplace_ellsort_double_real(
    void *scale_factor,
    bml_matrix_ellsort_t * A);

void bml_scale_inplace_ellsort_single_complex(
    void *scale_factor,
    bml_matrix_ellsort_t * A);

void bml_scale_inplace_ellsort_double_complex(
    void *scale_factor,
    bml_matrix_ellsort_t * A);

#endif
