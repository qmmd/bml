#ifndef __BML_DIAGONALIZE_DENSE_H
#define __BML_DIAGONALIZE_DENSE_H

#include "bml_types_dense.h"

void bml_diagonalize_dense(
    bml_matrix_dense_t * A,
    void *eigenvalues,
    bml_matrix_t * eigenvectors);

void bml_diagonalize_dense_single_real(
    bml_matrix_dense_t * A,
    void *eigenvalues,
    bml_matrix_dense_t * eigenvectors);

void bml_diagonalize_dense_double_real(
    bml_matrix_dense_t * A,
    void *eigenvalues,
    bml_matrix_dense_t * eigenvectors);

void bml_diagonalize_dense_single_complex(
    bml_matrix_dense_t * A,
    void *eigenvalues,
    bml_matrix_dense_t * eigenvectors);

void bml_diagonalize_dense_double_complex(
    bml_matrix_dense_t * A,
    void *eigenvalues,
    bml_matrix_dense_t * eigenvectors);

#endif
