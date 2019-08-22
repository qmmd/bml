/** \file */

#ifndef __BML_GETTERS_DENSE_H
#define __BML_GETTERS_DENSE_H

#include "bml_types_dense.h"

#include <complex.h>

void *bml_get_dense(
    bml_matrix_dense_t * A,
    int i,
    int j);

void *bml_get_dense_single_real(
    bml_matrix_dense_t * A,
    int i,
    int j);

void *bml_get_dense_double_real(
    bml_matrix_dense_t * A,
    int i,
    int j);

void *bml_get_dense_single_complex(
    bml_matrix_dense_t * A,
    int i,
    int j);

void *bml_get_dense_double_complex(
    bml_matrix_dense_t * A,
    int i,
    int j);

void *bml_get_row_dense(
    bml_matrix_dense_t * A,
    int i);

void *bml_get_row_dense_single_real(
    bml_matrix_dense_t * A,
    int i);

void *bml_get_row_dense_double_real(
    bml_matrix_dense_t * A,
    int i);

void *bml_get_row_dense_single_complex(
    bml_matrix_dense_t * A,
    int i);

void *bml_get_row_dense_double_complex(
    bml_matrix_dense_t * A,
    int i);

void *bml_get_diagonal_dense(
    bml_matrix_dense_t * A);

void *bml_get_diagonal_dense_single_real(
    bml_matrix_dense_t * A);

void *bml_get_diagonal_dense_double_real(
    bml_matrix_dense_t * A);

void *bml_get_diagonal_dense_single_complex(
    bml_matrix_dense_t * A);

void *bml_get_diagonal_dense_double_complex(
    bml_matrix_dense_t * A);

#endif
