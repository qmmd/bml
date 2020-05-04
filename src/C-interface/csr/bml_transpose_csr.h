#ifndef __BML_TRANSPOSE_CSR_H
#define __BML_TRANSPOSE_CSR_H

#include "bml_types_csr.h"

void csr_swap_row_entries_single_real(
    csr_sparse_row_t * row,
    const int ipos,
    const int jpos);

void csr_swap_row_entries_double_real(
    csr_sparse_row_t * row,
    const int ipos,
    const int jpos);

void csr_swap_row_entries_single_complex(
    csr_sparse_row_t * row,
    const int ipos,
    const int jpos);

void csr_swap_row_entries_double_complex(
    csr_sparse_row_t * row,
    const int ipos,
    const int jpos);

bml_matrix_csr_t *bml_transpose_new_csr(
    bml_matrix_csr_t * A);

bml_matrix_csr_t *bml_transpose_new_csr_single_real(
    bml_matrix_csr_t * A);

bml_matrix_csr_t *bml_transpose_new_csr_double_real(
    bml_matrix_csr_t * A);

bml_matrix_csr_t *bml_transpose_new_csr_single_complex(
    bml_matrix_csr_t * A);

bml_matrix_csr_t *bml_transpose_new_csr_double_complex(
    bml_matrix_csr_t * A);

void bml_transpose_csr(
    bml_matrix_csr_t * A);

void bml_transpose_csr_single_real(
    bml_matrix_csr_t * A);

void bml_transpose_csr_double_real(
    bml_matrix_csr_t * A);

void bml_transpose_csr_single_complex(
    bml_matrix_csr_t * A);

void bml_transpose_csr_double_complex(
    bml_matrix_csr_t * A);

#endif
