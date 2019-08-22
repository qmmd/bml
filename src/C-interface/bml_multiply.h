/** \file */

#ifndef __BML_MULTIPLY_H
#define __BML_MULTIPLY_H

#include "bml_types.h"

// Multiply - C = alpha * A * B + beta * C
void bml_multiply(
    bml_matrix_t * A,
    bml_matrix_t * B,
    bml_matrix_t * C,
    double alpha,
    double beta,
    double threshold);

// Multiply X^2 - X2 = X * X
void *bml_multiply_x2(
    bml_matrix_t * X,
    bml_matrix_t * X2,
    double threshold);

// Multiply - C = A * B
void bml_multiply_AB(
    bml_matrix_t * A,
    bml_matrix_t * B,
    bml_matrix_t * C,
    double threshold);

// Multiply with threshold adjustment - C = A * B
void bml_multiply_adjust_AB(
    bml_matrix_t * A,
    bml_matrix_t * B,
    bml_matrix_t * C,
    double threshold);

#endif
