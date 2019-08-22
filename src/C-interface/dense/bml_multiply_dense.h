#ifndef __BML_MULTIPLY_DENSE_H
#define __BML_MULTIPLY_DENSE_H

#include "bml_types_dense.h"

// Matrix multiply - C = alpha * A * B + beta * C
void bml_multiply_dense(
    bml_matrix_dense_t * A,
    bml_matrix_dense_t * B,
    bml_matrix_dense_t * C,
    double alpha,
    double beta);

void bml_multiply_dense_single_real(
    bml_matrix_dense_t * A,
    bml_matrix_dense_t * B,
    bml_matrix_dense_t * C,
    double alpha,
    double beta);

void bml_multiply_dense_double_real(
    bml_matrix_dense_t * A,
    bml_matrix_dense_t * B,
    bml_matrix_dense_t * C,
    double alpha,
    double beta);

void bml_multiply_dense_single_complex(
    bml_matrix_dense_t * A,
    bml_matrix_dense_t * B,
    bml_matrix_dense_t * C,
    double alpha,
    double beta);

void bml_multiply_dense_double_complex(
    bml_matrix_dense_t * A,
    bml_matrix_dense_t * B,
    bml_matrix_dense_t * C,
    double alpha,
    double beta);

void *bml_multiply_x2_dense(
    bml_matrix_dense_t * X,
    bml_matrix_dense_t * X2);

void *bml_multiply_x2_dense_single_real(
    bml_matrix_dense_t * X,
    bml_matrix_dense_t * X2);

void *bml_multiply_x2_dense_double_real(
    bml_matrix_dense_t * X,
    bml_matrix_dense_t * X2);

void *bml_multiply_x2_dense_single_complex(
    bml_matrix_dense_t * X,
    bml_matrix_dense_t * X2);

void *bml_multiply_x2_dense_double_complex(
    bml_matrix_dense_t * X,
    bml_matrix_dense_t * X2);

void bml_multiply_AB_dense(
    bml_matrix_dense_t * A,
    bml_matrix_dense_t * B,
    bml_matrix_dense_t * C);

void bml_multiply_AB_dense_single_real(
    bml_matrix_dense_t * A,
    bml_matrix_dense_t * B,
    bml_matrix_dense_t * C);

void bml_multiply_AB_dense_double_real(
    bml_matrix_dense_t * A,
    bml_matrix_dense_t * B,
    bml_matrix_dense_t * C);

void bml_multiply_AB_dense_single_complex(
    bml_matrix_dense_t * A,
    bml_matrix_dense_t * B,
    bml_matrix_dense_t * C);

void bml_multiply_AB_dense_double_complex(
    bml_matrix_dense_t * A,
    bml_matrix_dense_t * B,
    bml_matrix_dense_t * C);

void bml_multiply_adjust_AB_dense(
    bml_matrix_dense_t * A,
    bml_matrix_dense_t * B,
    bml_matrix_dense_t * C);

void bml_multiply_adjust_AB_dense_single_real(
    bml_matrix_dense_t * A,
    bml_matrix_dense_t * B,
    bml_matrix_dense_t * C);

void bml_multiply_adjust_AB_dense_double_real(
    bml_matrix_dense_t * A,
    bml_matrix_dense_t * B,
    bml_matrix_dense_t * C);

void bml_multiply_adjust_AB_dense_single_complex(
    bml_matrix_dense_t * A,
    bml_matrix_dense_t * B,
    bml_matrix_dense_t * C);

void bml_multiply_adjust_AB_dense_double_complex(
    bml_matrix_dense_t * A,
    bml_matrix_dense_t * B,
    bml_matrix_dense_t * C);

#endif
