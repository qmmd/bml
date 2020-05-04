#ifndef __BML_DIAGONALIZE_CSR_H
#define __BML_DIAGONALIZE_CSR_H

#include "bml_types_csr.h"

void bml_diagonalize_csr(
    bml_matrix_csr_t * A,
    void *eigenvalues,
    bml_matrix_t * eigenvectors);

void bml_diagonalize_csr_single_real(
    bml_matrix_csr_t * A,
    void *eigenvalues,
    bml_matrix_csr_t * eigenvectors);

void bml_diagonalize_csr_double_real(
    bml_matrix_csr_t * A,
    void *eigenvalues,
    bml_matrix_csr_t * eigenvectors);

void bml_diagonalize_csr_single_complex(
    bml_matrix_csr_t * A,
    void *eigenvalues,
    bml_matrix_csr_t * eigenvectors);

void bml_diagonalize_csr_double_complex(
    bml_matrix_csr_t * A,
    void *eigenvalues,
    bml_matrix_csr_t * eigenvectors);

#endif
