#ifndef __BML_NORMALIZE_ELLSORT_H
#define __BML_NORMALIZE_ELLSORT_H

#include "bml_types_ellsort.h"

void bml_normalize_ellsort(
    bml_matrix_ellsort_t * A,
    double mineval,
    double maxeval);

void bml_normalize_ellsort_single_real(
    bml_matrix_ellsort_t * A,
    double mineval,
    double maxeval);

void bml_normalize_ellsort_double_real(
    bml_matrix_ellsort_t * A,
    double mineval,
    double maxeval);

void bml_normalize_ellsort_single_complex(
    bml_matrix_ellsort_t * A,
    double mineval,
    double maxeval);

void bml_normalize_ellsort_double_complex(
    bml_matrix_ellsort_t * A,
    double mineval,
    double maxeval);

void *bml_gershgorin_ellsort(
    bml_matrix_ellsort_t * A);

void *bml_gershgorin_ellsort_single_real(
    bml_matrix_ellsort_t * A);

void *bml_gershgorin_ellsort_double_real(
    bml_matrix_ellsort_t * A);

void *bml_gershgorin_ellsort_single_complex(
    bml_matrix_ellsort_t * A);

void *bml_gershgorin_ellsort_double_complex(
    bml_matrix_ellsort_t * A);

void *bml_gershgorin_partial_ellsort(
    bml_matrix_ellsort_t * A,
    int nrows);

void *bml_gershgorin_partial_ellsort_single_real(
    bml_matrix_ellsort_t * A,
    int nrows);

void *bml_gershgorin_partial_ellsort_double_real(
    bml_matrix_ellsort_t * A,
    int nrows);

void *bml_gershgorin_partial_ellsort_single_complex(
    bml_matrix_ellsort_t * A,
    int nrows);

void *bml_gershgorin_partial_ellsort_double_complex(
    bml_matrix_ellsort_t * A,
    int nrows);

#endif
