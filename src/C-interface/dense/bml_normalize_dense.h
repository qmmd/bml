#ifndef __BML_NORMALIZE_DENSE_H
#define __BML_NORMALIZE_DENSE_H

#include "bml_types_dense.h"

void bml_normalize_dense(
    bml_matrix_dense_t * A,
    double mineval,
    double maxeval);

void bml_normalize_dense_single_real(
    bml_matrix_dense_t * A,
    double maxeval,
    double maxminusmin);

void bml_normalize_dense_double_real(
    bml_matrix_dense_t * A,
    double mineval,
    double maxeval);

void bml_normalize_dense_single_complex(
    bml_matrix_dense_t * A,
    double mineval,
    double maxeval);

void bml_normalize_dense_double_complex(
    bml_matrix_dense_t * A,
    double mineval,
    double maxeval);

void *bml_gershgorin_dense(
    bml_matrix_dense_t * A);

void *bml_gershgorin_dense_single_real(
    bml_matrix_dense_t * A);

void *bml_gershgorin_dense_double_real(
    bml_matrix_dense_t * A);

void *bml_gershgorin_dense_single_complex(
    bml_matrix_dense_t * A);

void *bml_gershgorin_dense_double_complex(
    bml_matrix_dense_t * A);

void *bml_gershgorin_partial_dense(
    bml_matrix_dense_t * A,
    int nrows);

void *bml_gershgorin_partial_dense_single_real(
    bml_matrix_dense_t * A,
    int nrows);

void *bml_gershgorin_partial_dense_double_real(
    bml_matrix_dense_t * A,
    int nrows);

void *bml_gershgorin_partial_dense_single_complex(
    bml_matrix_dense_t * A,
    int nrows);

void *bml_gershgorin_partial_dense_double_complex(
    bml_matrix_dense_t * A,
    int nrows);

#endif
