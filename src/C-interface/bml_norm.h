/** \file */

#ifndef __BML_NORM_H
#define __BML_NORM_H

#include "bml_types.h"

// Calculate the sum of squares of all the elements in A
double bml_sum_squares(
    bml_matrix_t * A);

// Calculate the sum of squares of all the elements of
// alpha * A + beta * B
double bml_sum_squares2(
    bml_matrix_t * A,
    bml_matrix_t * B,
    double alpha,
    double beta,
    double threshold);

// Calculate the sum of squares for submatrix core elements
double bml_sum_squares_submatrix(
    bml_matrix_t * A,
    int core_size);

// Calculate Frobenius norm for matrix A
// sqrt(sum(A_ij*A_ij)
double bml_fnorm(
    bml_matrix_t * A);

// Calculate Frobenius norm for 2 matrices
double bml_fnorm2(
    bml_matrix_t * A,
    bml_matrix_t * B);

#endif
