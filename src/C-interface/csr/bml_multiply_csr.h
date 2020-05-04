#ifndef __BML_MULTIPLY_CSR_H
#define __BML_MULTIPLY_CSR_H

#include "bml_types_csr.h"

void bml_multiply_csr(
    bml_matrix_csr_t * A,
    bml_matrix_csr_t * B,
    bml_matrix_csr_t * C,
    double alpha,
    double beta,
    double threshold);

void bml_multiply_csr_single_real(
    bml_matrix_csr_t * A,
    bml_matrix_csr_t * B,
    bml_matrix_csr_t * C,
    double alpha,
    double beta,
    double threshold);

void bml_multiply_csr_double_real(
    bml_matrix_csr_t * A,
    bml_matrix_csr_t * B,
    bml_matrix_csr_t * C,
    double alpha,
    double beta,
    double threshold);

void bml_multiply_csr_single_complex(
    bml_matrix_csr_t * A,
    bml_matrix_csr_t * B,
    bml_matrix_csr_t * C,
    double alpha,
    double beta,
    double threshold);

void bml_multiply_csr_double_complex(
    bml_matrix_csr_t * A,
    bml_matrix_csr_t * B,
    bml_matrix_csr_t * C,
    double alpha,
    double beta,
    double threshold);

void *bml_multiply_x2_csr(
    bml_matrix_csr_t * X,
    bml_matrix_csr_t * X2,
    double threshold);

void *bml_multiply_x2_csr_single_real(
    bml_matrix_csr_t * X,
    bml_matrix_csr_t * X2,
    double threshold);

void *bml_multiply_x2_csr_double_real(
    bml_matrix_csr_t * X,
    bml_matrix_csr_t * X2,
    double threshold);

void *bml_multiply_x2_csr_single_complex(
    bml_matrix_csr_t * X,
    bml_matrix_csr_t * X2,
    double threshold);

void *bml_multiply_x2_csr_double_complex(
    bml_matrix_csr_t * X,
    bml_matrix_csr_t * X2,
    double threshold);

void bml_multiply_AB_csr(
    bml_matrix_csr_t * A,
    bml_matrix_csr_t * B,
    bml_matrix_csr_t * C,
    double threshold);

void bml_multiply_AB_csr_single_real(
    bml_matrix_csr_t * A,
    bml_matrix_csr_t * B,
    bml_matrix_csr_t * C,
    double threshold);

void bml_multiply_AB_csr_double_real(
    bml_matrix_csr_t * A,
    bml_matrix_csr_t * B,
    bml_matrix_csr_t * C,
    double threshold);

void bml_multiply_AB_csr_single_complex(
    bml_matrix_csr_t * A,
    bml_matrix_csr_t * B,
    bml_matrix_csr_t * C,
    double threshold);

void bml_multiply_AB_csr_double_complex(
    bml_matrix_csr_t * A,
    bml_matrix_csr_t * B,
    bml_matrix_csr_t * C,
    double threshold);

void bml_multiply_adjust_AB_csr(
    bml_matrix_csr_t * A,
    bml_matrix_csr_t * B,
    bml_matrix_csr_t * C,
    double threshold);

void bml_multiply_adjust_AB_csr_single_real(
    bml_matrix_csr_t * A,
    bml_matrix_csr_t * B,
    bml_matrix_csr_t * C,
    double threshold);

void bml_multiply_adjust_AB_csr_double_real(
    bml_matrix_csr_t * A,
    bml_matrix_csr_t * B,
    bml_matrix_csr_t * C,
    double threshold);

void bml_multiply_adjust_AB_csr_single_complex(
    bml_matrix_csr_t * A,
    bml_matrix_csr_t * B,
    bml_matrix_csr_t * C,
    double threshold);

void bml_multiply_adjust_AB_csr_double_complex(
    bml_matrix_csr_t * A,
    bml_matrix_csr_t * B,
    bml_matrix_csr_t * C,
    double threshold);

#endif
