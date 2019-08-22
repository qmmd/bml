#ifndef __BML_NORM_DENSE_H
#define __BML_NORM_DENSE_H

#include "bml_types_dense.h"

double bml_sum_squares_dense(
    bml_matrix_dense_t * A);

double bml_sum_squares_dense_single_real(
    bml_matrix_dense_t * A);

double bml_sum_squares_dense_double_real(
    bml_matrix_dense_t * A);

double bml_sum_squares_dense_single_complex(
    bml_matrix_dense_t * A);

double bml_sum_squares_dense_double_complex(
    bml_matrix_dense_t * A);

double bml_sum_squares_submatrix_dense(
    bml_matrix_dense_t * A,
    int core_size);

double bml_sum_squares_submatrix_dense_single_real(
    bml_matrix_dense_t * A,
    int core_size);

double bml_sum_squares_submatrix_dense_double_real(
    bml_matrix_dense_t * A,
    int core_size);

double bml_sum_squares_submatrix_dense_single_complex(
    bml_matrix_dense_t * A,
    int core_size);

double bml_sum_squares_submatrix_dense_double_complex(
    bml_matrix_dense_t * A,
    int core_size);

double bml_sum_squares2_dense(
    bml_matrix_dense_t * A,
    bml_matrix_dense_t * B,
    double alpha,
    double beta,
    double threshold);

double bml_sum_squares2_dense_single_real(
    bml_matrix_dense_t * A,
    bml_matrix_dense_t * B,
    double alpha,
    double beta,
    double threshold);

double bml_sum_squares2_dense_double_real(
    bml_matrix_dense_t * A,
    bml_matrix_dense_t * B,
    double alpha,
    double beta,
    double threshold);

double bml_sum_squares2_dense_single_complex(
    bml_matrix_dense_t * A,
    bml_matrix_dense_t * B,
    double alpha,
    double beta,
    double threshold);

double bml_sum_squares2_dense_double_complex(
    bml_matrix_dense_t * A,
    bml_matrix_dense_t * B,
    double alpha,
    double beta,
    double threshold);

double bml_fnorm_dense(
    bml_matrix_dense_t * A);

double bml_fnorm_dense_single_real(
    bml_matrix_dense_t * A);

double bml_fnorm_dense_double_real(
    bml_matrix_dense_t * A);

double bml_fnorm_dense_single_complex(
    bml_matrix_dense_t * A);

double bml_fnorm_dense_double_complex(
    bml_matrix_dense_t * A);

double bml_fnorm2_dense(
    bml_matrix_dense_t * A,
    bml_matrix_dense_t * B);

double bml_fnorm2_dense_single_real(
    bml_matrix_dense_t * A,
    bml_matrix_dense_t * B);

double bml_fnorm2_dense_double_real(
    bml_matrix_dense_t * A,
    bml_matrix_dense_t * B);

double bml_fnorm2_dense_single_complex(
    bml_matrix_dense_t * A,
    bml_matrix_dense_t * B);

double bml_fnorm2_dense_double_complex(
    bml_matrix_dense_t * A,
    bml_matrix_dense_t * B);

#endif
