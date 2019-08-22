#ifndef __BML_NORM_ELLSORT_H
#define __BML_NORM_ELLSORT_H

#include "bml_types_ellsort.h"

double bml_sum_squares_ellsort(
    bml_matrix_ellsort_t * A);

double bml_sum_squares_ellsort_single_real(
    bml_matrix_ellsort_t * A);

double bml_sum_squares_ellsort_double_real(
    bml_matrix_ellsort_t * A);

double bml_sum_squares_ellsort_single_complex(
    bml_matrix_ellsort_t * A);

double bml_sum_squares_ellsort_double_complex(
    bml_matrix_ellsort_t * A);

double bml_sum_squares_submatrix_ellsort(
    bml_matrix_ellsort_t * A,
    int core_size);

double bml_sum_squares_submatrix_ellsort_single_real(
    bml_matrix_ellsort_t * A,
    int core_size);

double bml_sum_squares_submatrix_ellsort_double_real(
    bml_matrix_ellsort_t * A,
    int core_size);

double bml_sum_squares_submatrix_ellsort_single_complex(
    bml_matrix_ellsort_t * A,
    int core_size);

double bml_sum_squares_submatrix_ellsort_double_complex(
    bml_matrix_ellsort_t * A,
    int core_size);

double bml_sum_squares2_ellsort(
    bml_matrix_ellsort_t * A,
    bml_matrix_ellsort_t * B,
    double alpha,
    double beta,
    double threshold);

double bml_sum_squares2_ellsort_single_real(
    bml_matrix_ellsort_t * A,
    bml_matrix_ellsort_t * B,
    double alpha,
    double beta,
    double threshold);

double bml_sum_squares2_ellsort_double_real(
    bml_matrix_ellsort_t * A,
    bml_matrix_ellsort_t * B,
    double alpha,
    double beta,
    double threshold);

double bml_sum_squares2_ellsort_single_complex(
    bml_matrix_ellsort_t * A,
    bml_matrix_ellsort_t * B,
    double alpha,
    double beta,
    double threshold);

double bml_sum_squares2_ellsort_double_complex(
    bml_matrix_ellsort_t * A,
    bml_matrix_ellsort_t * B,
    double alpha,
    double beta,
    double threshold);

double bml_fnorm_ellsort(
    bml_matrix_ellsort_t * A);

double bml_fnorm_ellsort_single_real(
    bml_matrix_ellsort_t * A);

double bml_fnorm_ellsort_double_real(
    bml_matrix_ellsort_t * A);

double bml_fnorm_ellsort_single_complex(
    bml_matrix_ellsort_t * A);

double bml_fnorm_ellsort_double_complex(
    bml_matrix_ellsort_t * A);

double bml_fnorm2_ellsort(
    bml_matrix_ellsort_t * A,
    bml_matrix_ellsort_t * B);

double bml_fnorm2_ellsort_single_real(
    bml_matrix_ellsort_t * A,
    bml_matrix_ellsort_t * B);

double bml_fnorm2_ellsort_double_real(
    bml_matrix_ellsort_t * A,
    bml_matrix_ellsort_t * B);

double bml_fnorm2_ellsort_single_complex(
    bml_matrix_ellsort_t * A,
    bml_matrix_ellsort_t * B);

double bml_fnorm2_ellsort_double_complex(
    bml_matrix_ellsort_t * A,
    bml_matrix_ellsort_t * B);

#endif
