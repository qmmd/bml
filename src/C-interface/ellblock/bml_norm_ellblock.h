#ifndef __BML_NORM_ELLBLOCK_H
#define __BML_NORM_ELLBLOCK_H

#include "bml_types_ellblock.h"

double bml_sum_squares_ellblock(
    bml_matrix_ellblock_t * A);

double bml_sum_squares_ellblock_single_real(
    bml_matrix_ellblock_t * A);

double bml_sum_squares_ellblock_double_real(
    bml_matrix_ellblock_t * A);

double bml_sum_squares_ellblock_single_complex(
    bml_matrix_ellblock_t * A);

double bml_sum_squares_ellblock_double_complex(
    bml_matrix_ellblock_t * A);

double bml_sum_squares_submatrix_ellblock(
    bml_matrix_ellblock_t * A,
    int core_size);

double bml_sum_squares_submatrix_ellblock_single_real(
    bml_matrix_ellblock_t * A,
    int core_size);

double bml_sum_squares_submatrix_ellblock_double_real(
    bml_matrix_ellblock_t * A,
    int core_size);

double bml_sum_squares_submatrix_ellblock_single_complex(
    bml_matrix_ellblock_t * A,
    int core_size);

double bml_sum_squares_submatrix_ellblock_double_complex(
    bml_matrix_ellblock_t * A,
    int core_size);

double bml_sum_squares2_ellblock(
    bml_matrix_ellblock_t * A,
    bml_matrix_ellblock_t * B,
    double alpha,
    double beta,
    double threshold);

double bml_sum_squares2_ellblock_single_real(
    bml_matrix_ellblock_t * A,
    bml_matrix_ellblock_t * B,
    double alpha,
    double beta,
    double threshold);

double bml_sum_squares2_ellblock_double_real(
    bml_matrix_ellblock_t * A,
    bml_matrix_ellblock_t * B,
    double alpha,
    double beta,
    double threshold);

double bml_sum_squares2_ellblock_single_complex(
    bml_matrix_ellblock_t * A,
    bml_matrix_ellblock_t * B,
    double alpha,
    double beta,
    double threshold);

double bml_sum_squares2_ellblock_double_complex(
    bml_matrix_ellblock_t * A,
    bml_matrix_ellblock_t * B,
    double alpha,
    double beta,
    double threshold);

double bml_fnorm_ellblock(
    bml_matrix_ellblock_t * A);

double bml_fnorm_ellblock_single_real(
    bml_matrix_ellblock_t * A);

double bml_fnorm_ellblock_double_real(
    bml_matrix_ellblock_t * A);

double bml_fnorm_ellblock_single_complex(
    bml_matrix_ellblock_t * A);

double bml_fnorm_ellblock_double_complex(
    bml_matrix_ellblock_t * A);

double bml_fnorm2_ellblock(
    bml_matrix_ellblock_t * A,
    bml_matrix_ellblock_t * B);

double bml_fnorm2_ellblock_single_real(
    bml_matrix_ellblock_t * A,
    bml_matrix_ellblock_t * B);

double bml_fnorm2_ellblock_double_real(
    bml_matrix_ellblock_t * A,
    bml_matrix_ellblock_t * B);

double bml_fnorm2_ellblock_single_complex(
    bml_matrix_ellblock_t * A,
    bml_matrix_ellblock_t * B);

double bml_fnorm2_ellblock_double_complex(
    bml_matrix_ellblock_t * A,
    bml_matrix_ellblock_t * B);

#endif
