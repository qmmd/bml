/** \file */

#ifndef __BML_SETTERS_CSR_H
#define __BML_SETTERS_CSR_H

#include "bml_types_csr.h"

#include <complex.h>
/*
void csr_set_row_element_new(
    csr_sparse_row_t * arow,
    const int j,
    const void *element);
*/
void csr_set_row_element_new_single_real(
    csr_sparse_row_t * arow,
    const int j,
    const void *element);

void csr_set_row_element_new_double_real(
    csr_sparse_row_t * arow,
    const int j,
    const void *element);

void csr_set_row_element_new_single_complex(
    csr_sparse_row_t * arow,
    const int j,
    const void *element);

void csr_set_row_element_new_double_complex(
    csr_sparse_row_t * arow,
    const int j,
    const void *element);
/*
void csr_set_row_element(
    csr_sparse_row_t * arow,
    const int j,
    const void *element);
*/

void csr_set_row_element_single_real(
    csr_sparse_row_t * arow,
    const int j,
    const void *element);

void csr_set_row_element_double_real(
    csr_sparse_row_t * arow,
    const int j,
    const void *element);

void csr_set_row_element_single_complex(
    csr_sparse_row_t * arow,
    const int j,
    const void *element);

void csr_set_row_element_double_complex(
    csr_sparse_row_t * arow,
    const int j,
    const void *element);
/*
void csr_set_row(
    csr_sparse_row_t * arow,
    const int count,
    const int * cols,
    const void * vals,
    const double threshold);
*/

void csr_set_row_single_real(
    csr_sparse_row_t * arow,
    const int count,
    void *vals,
    const double threshold);

void csr_set_row_double_real(
    csr_sparse_row_t * arow,
    const int count,
    void *vals,
    const double threshold);

void csr_set_row_single_complex(
    csr_sparse_row_t * arow,
    const int count,
    void *vals,
    const double threshold);

void csr_set_row_double_complex(
    csr_sparse_row_t * arow,
    const int count,
    void *vals,
    const double threshold);

void csr_set_sparse_row_single_real(
    csr_sparse_row_t * arow,
    const int count,
    const int *cols,
    void *vals,
    const double threshold);

void csr_set_sparse_row_double_real(
    csr_sparse_row_t * arow,
    const int count,
    const int *cols,
    void *vals,
    const double threshold);

void csr_set_sparse_row_single_complex(
    csr_sparse_row_t * arow,
    const int count,
    const int *cols,
    void *vals,
    const double threshold);

void csr_set_sparse_row_double_complex(
    csr_sparse_row_t * arow,
    const int count,
    const int *cols,
    void *vals,
    const double threshold);

void bml_set_element_new_csr(
    bml_matrix_csr_t * A,
    const int i,
    const int j,
    void *value);

void bml_set_element_new_csr_single_real(
    bml_matrix_csr_t * A,
    const int i,
    const int j,
    void *value);

void bml_set_element_new_csr_double_real(
    bml_matrix_csr_t * A,
    const int i,
    const int j,
    void *value);

void bml_set_element_new_csr_single_complex(
    bml_matrix_csr_t * A,
    const int i,
    const int j,
    void *value);

void bml_set_element_new_csr_double_complex(
    bml_matrix_csr_t * A,
    const int i,
    const int j,
    void *value);

void bml_set_element_csr(
    bml_matrix_csr_t * A,
    const int i,
    const int j,
    void *value);

void bml_set_element_csr_single_real(
    bml_matrix_csr_t * A,
    const int i,
    const int j,
    void *value);

void bml_set_element_csr_double_real(
    bml_matrix_csr_t * A,
    const int i,
    const int j,
    void *value);

void bml_set_element_csr_single_complex(
    bml_matrix_csr_t * A,
    const int i,
    const int j,
    void *value);

void bml_set_element_csr_double_complex(
    bml_matrix_csr_t * A,
    const int i,
    const int j,
    void *value);

void bml_set_row_csr(
    bml_matrix_csr_t * A,
    const int i,
    void *row,
    const double threshold);

void bml_set_row_csr_single_real(
    bml_matrix_csr_t * A,
    const int i,
    void *row,
    const double threshold);

void bml_set_row_csr_double_real(
    bml_matrix_csr_t * A,
    const int i,
    void *row,
    const double threshold);

void bml_set_row_csr_single_complex(
    bml_matrix_csr_t * A,
    const int i,
    void *row,
    const double threshold);

void bml_set_row_csr_double_complex(
    bml_matrix_csr_t * A,
    const int i,
    void *row,
    const double threshold);

void bml_set_diagonal_csr(
    bml_matrix_csr_t * A,
    void *diagonal,
    const double threshold);

void bml_set_diagonal_csr_single_real(
    bml_matrix_csr_t * A,
    void *diagonal,
    const double threshold);

void bml_set_diagonal_csr_double_real(
    bml_matrix_csr_t * A,
    void *diagonal,
    const double threshold);

void bml_set_diagonal_csr_single_complex(
    bml_matrix_csr_t * A,
    void *diagonal,
    const double threshold);

void bml_set_diagonal_csr_double_complex(
    bml_matrix_csr_t * A,
    void *diagonal,
    const double threshold);


void bml_set_sparse_row_csr(
    bml_matrix_csr_t * A,
    const int i,
    const int count,
    const int *cols,
    void *vals,
    const double threshold);

void bml_set_sparse_row_csr_single_real(
    bml_matrix_csr_t * A,
    const int i,
    const int count,
    const int *cols,
    void *vals,
    const double threshold);

void bml_set_sparse_row_csr_double_real(
    bml_matrix_csr_t * A,
    const int i,
    const int count,
    const int *cols,
    void *vals,
    const double threshold);

void bml_set_sparse_row_csr_single_complex(
    bml_matrix_csr_t * A,
    const int i,
    const int count,
    const int *cols,
    void *vals,
    const double threshold);

void bml_set_sparse_row_csr_double_complex(
    bml_matrix_csr_t * A,
    const int i,
    const int count,
    const int *cols,
    void *vals,
    const double threshold);

#endif
