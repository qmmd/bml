#ifndef __BML_NORM_CSR_H
#define __BML_NORM_CSR_H

#include "bml_types_csr.h"

double bml_sum_squares_csr(
    bml_matrix_csr_t * A);

double bml_sum_squares_csr_single_real(
    bml_matrix_csr_t * A);

double bml_sum_squares_csr_double_real(
    bml_matrix_csr_t * A);

double bml_sum_squares_csr_single_complex(
    bml_matrix_csr_t * A);

double bml_sum_squares_csr_double_complex(
    bml_matrix_csr_t * A);

double bml_sum_squares_submatrix_csr(
    bml_matrix_csr_t * A,
    int core_size);

double bml_sum_squares_submatrix_csr_single_real(
    bml_matrix_csr_t * A,
    int core_size);

double bml_sum_squares_submatrix_csr_double_real(
    bml_matrix_csr_t * A,
    int core_size);

double bml_sum_squares_submatrix_csr_single_complex(
    bml_matrix_csr_t * A,
    int core_size);

double bml_sum_squares_submatrix_csr_double_complex(
    bml_matrix_csr_t * A,
    int core_size);

double bml_sum_squares2_csr(
    bml_matrix_csr_t * A,
    bml_matrix_csr_t * B,
    double alpha,
    double beta,
    double threshold);

double bml_sum_squares2_csr_single_real(
    bml_matrix_csr_t * A,
    bml_matrix_csr_t * B,
    double alpha,
    double beta,
    double threshold);

double bml_sum_squares2_csr_double_real(
    bml_matrix_csr_t * A,
    bml_matrix_csr_t * B,
    double alpha,
    double beta,
    double threshold);

double bml_sum_squares2_csr_single_complex(
    bml_matrix_csr_t * A,
    bml_matrix_csr_t * B,
    double alpha,
    double beta,
    double threshold);

double bml_sum_squares2_csr_double_complex(
    bml_matrix_csr_t * A,
    bml_matrix_csr_t * B,
    double alpha,
    double beta,
    double threshold);

double bml_fnorm_csr(
    bml_matrix_csr_t * A);

double bml_fnorm_csr_single_real(
    bml_matrix_csr_t * A);

double bml_fnorm_csr_double_real(
    bml_matrix_csr_t * A);

double bml_fnorm_csr_single_complex(
    bml_matrix_csr_t * A);

double bml_fnorm_csr_double_complex(
    bml_matrix_csr_t * A);

double bml_fnorm2_csr(
    bml_matrix_csr_t * A,
    bml_matrix_csr_t * B);

double bml_fnorm2_csr_single_real(
    bml_matrix_csr_t * A,
    bml_matrix_csr_t * B);

double bml_fnorm2_csr_double_real(
    bml_matrix_csr_t * A,
    bml_matrix_csr_t * B);

double bml_fnorm2_csr_single_complex(
    bml_matrix_csr_t * A,
    bml_matrix_csr_t * B);

double bml_fnorm2_csr_double_complex(
    bml_matrix_csr_t * A,
    bml_matrix_csr_t * B);

#endif
