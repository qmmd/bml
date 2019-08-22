#ifndef __BML_NORM_ELLPACK_H
#define __BML_NORM_ELLPACK_H

#include "bml_types_ellpack.h"

double bml_sum_squares_ellpack(
    bml_matrix_ellpack_t * A);

double bml_sum_squares_ellpack_single_real(
    bml_matrix_ellpack_t * A);

double bml_sum_squares_ellpack_double_real(
    bml_matrix_ellpack_t * A);

double bml_sum_squares_ellpack_single_complex(
    bml_matrix_ellpack_t * A);

double bml_sum_squares_ellpack_double_complex(
    bml_matrix_ellpack_t * A);

double bml_sum_squares_submatrix_ellpack(
    bml_matrix_ellpack_t * A,
    int core_size);

double bml_sum_squares_submatrix_ellpack_single_real(
    bml_matrix_ellpack_t * A,
    int core_size);

double bml_sum_squares_submatrix_ellpack_double_real(
    bml_matrix_ellpack_t * A,
    int core_size);

double bml_sum_squares_submatrix_ellpack_single_complex(
    bml_matrix_ellpack_t * A,
    int core_size);

double bml_sum_squares_submatrix_ellpack_double_complex(
    bml_matrix_ellpack_t * A,
    int core_size);

double bml_sum_squares2_ellpack(
    bml_matrix_ellpack_t * A,
    bml_matrix_ellpack_t * B,
    double alpha,
    double beta,
    double threshold);

double bml_sum_squares2_ellpack_single_real(
    bml_matrix_ellpack_t * A,
    bml_matrix_ellpack_t * B,
    double alpha,
    double beta,
    double threshold);

double bml_sum_squares2_ellpack_double_real(
    bml_matrix_ellpack_t * A,
    bml_matrix_ellpack_t * B,
    double alpha,
    double beta,
    double threshold);

double bml_sum_squares2_ellpack_single_complex(
    bml_matrix_ellpack_t * A,
    bml_matrix_ellpack_t * B,
    double alpha,
    double beta,
    double threshold);

double bml_sum_squares2_ellpack_double_complex(
    bml_matrix_ellpack_t * A,
    bml_matrix_ellpack_t * B,
    double alpha,
    double beta,
    double threshold);

double bml_fnorm_ellpack(
    bml_matrix_ellpack_t * A);

double bml_fnorm_ellpack_single_real(
    bml_matrix_ellpack_t * A);

double bml_fnorm_ellpack_double_real(
    bml_matrix_ellpack_t * A);

double bml_fnorm_ellpack_single_complex(
    bml_matrix_ellpack_t * A);

double bml_fnorm_ellpack_double_complex(
    bml_matrix_ellpack_t * A);

double bml_fnorm2_ellpack(
    bml_matrix_ellpack_t * A,
    bml_matrix_ellpack_t * B);

double bml_fnorm2_ellpack_single_real(
    bml_matrix_ellpack_t * A,
    bml_matrix_ellpack_t * B);

double bml_fnorm2_ellpack_double_real(
    bml_matrix_ellpack_t * A,
    bml_matrix_ellpack_t * B);

double bml_fnorm2_ellpack_single_complex(
    bml_matrix_ellpack_t * A,
    bml_matrix_ellpack_t * B);

double bml_fnorm2_ellpack_double_complex(
    bml_matrix_ellpack_t * A,
    bml_matrix_ellpack_t * B);

#endif
