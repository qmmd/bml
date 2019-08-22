/** \file */

#ifndef __BML_GETTERS_ELLSORT_H
#define __BML_GETTERS_ELLSORT_H

#include "bml_types_ellsort.h"

#include <complex.h>

void *bml_get_ellsort(
    bml_matrix_ellsort_t * A,
    int i,
    int j);

void *bml_get_ellsort_single_real(
    bml_matrix_ellsort_t * A,
    int i,
    int j);

void *bml_get_ellsort_double_real(
    bml_matrix_ellsort_t * A,
    int i,
    int j);

void *bml_get_ellsort_single_complex(
    bml_matrix_ellsort_t * A,
    int i,
    int j);

void *bml_get_ellsort_double_complex(
    bml_matrix_ellsort_t * A,
    int i,
    int j);

void *bml_get_row_ellsort(
    bml_matrix_ellsort_t * A,
    int i);

void *bml_get_row_ellsort_single_real(
    bml_matrix_ellsort_t * A,
    int i);

void *bml_get_row_ellsort_double_real(
    bml_matrix_ellsort_t * A,
    int i);

void *bml_get_row_ellsort_single_complex(
    bml_matrix_ellsort_t * A,
    int i);

void *bml_get_row_ellsort_double_complex(
    bml_matrix_ellsort_t * A,
    int i);

void *bml_get_diagonal_ellsort(
    bml_matrix_ellsort_t * A);

void *bml_get_diagonal_ellsort_single_real(
    bml_matrix_ellsort_t * A);

void *bml_get_diagonal_ellsort_double_real(
    bml_matrix_ellsort_t * A);

void *bml_get_diagonal_ellsort_single_complex(
    bml_matrix_ellsort_t * A);

void *bml_get_diagonal_ellsort_double_complex(
    bml_matrix_ellsort_t * A);

#endif
