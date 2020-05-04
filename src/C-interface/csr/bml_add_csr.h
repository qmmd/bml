#ifndef __BML_ADD_CSR_H
#define __BML_ADD_CSR_H

#include "bml_types_csr.h"

void bml_add_csr(
    bml_matrix_csr_t * A,
    bml_matrix_csr_t * B,
    double alpha,
    double beta,
    double threshold);

void bml_add_csr_single_real(
    bml_matrix_csr_t * A,
    bml_matrix_csr_t * B,
    double alpha,
    double beta,
    double threshold);

void bml_add_csr_double_real(
    bml_matrix_csr_t * A,
    bml_matrix_csr_t * B,
    double alpha,
    double beta,
    double threshold);

void bml_add_csr_single_complex(
    bml_matrix_csr_t * A,
    bml_matrix_csr_t * B,
    double alpha,
    double beta,
    double threshold);

void bml_add_csr_double_complex(
    bml_matrix_csr_t * A,
    bml_matrix_csr_t * B,
    double alpha,
    double beta,
    double threshold);

double bml_add_norm_csr(
    bml_matrix_csr_t * A,
    bml_matrix_csr_t * B,
    double alpha,
    double beta,
    double threshold);

double bml_add_norm_csr_single_real(
    bml_matrix_csr_t * A,
    bml_matrix_csr_t * B,
    double alpha,
    double beta,
    double threshold);

double bml_add_norm_csr_double_real(
    bml_matrix_csr_t * A,
    bml_matrix_csr_t * B,
    double alpha,
    double beta,
    double threshold);

double bml_add_norm_csr_single_complex(
    bml_matrix_csr_t * A,
    bml_matrix_csr_t * B,
    double alpha,
    double beta,
    double threshold);

double bml_add_norm_csr_double_complex(
    bml_matrix_csr_t * A,
    bml_matrix_csr_t * B,
    double alpha,
    double beta,
    double threshold);

void bml_add_identity_csr(
    bml_matrix_csr_t * A,
    double beta,
    double threshold);

void bml_add_identity_csr_single_real(
    bml_matrix_csr_t * A,
    double beta,
    double threshold);

void bml_add_identity_csr_double_real(
    bml_matrix_csr_t * A,
    double beta,
    double threshold);

void bml_add_identity_csr_single_complex(
    bml_matrix_csr_t * A,
    double beta,
    double threshold);

void bml_add_identity_csr_double_complex(
    bml_matrix_csr_t * A,
    double beta,
    double threshold);

void bml_scale_add_identity_csr(
    bml_matrix_csr_t * A,
    double alpha,
    double beta,
    double threshold);

void bml_scale_add_identity_csr_single_real(
    bml_matrix_csr_t * A,
    double alpha,
    double beta,
    double threshold);

void bml_scale_add_identity_csr_double_real(
    bml_matrix_csr_t * A,
    double alpha,
    double beta,
    double threshold);

void bml_scale_add_identity_csr_single_complex(
    bml_matrix_csr_t * A,
    double alpha,
    double beta,
    double threshold);

void bml_scale_add_identity_csr_double_complex(
    bml_matrix_csr_t * A,
    double alpha,
    double beta,
    double threshold);

#endif
