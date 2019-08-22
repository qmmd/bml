#ifndef __BML_DIAGONALIZE_BELLBLOCK_H
#define __BML_DIAGONALIZE_BELLBLOCK_H

#include "bml_types_ellblock.h"

void bml_diagonalize_ellblock(
    bml_matrix_ellblock_t * A,
    void *eigenvalues,
    bml_matrix_t * eigenvectors);

void bml_diagonalize_ellblock_single_real(
    bml_matrix_ellblock_t * A,
    void *eigenvalues,
    bml_matrix_ellblock_t * eigenvectors);

void bml_diagonalize_ellblock_double_real(
    bml_matrix_ellblock_t * A,
    void *eigenvalues,
    bml_matrix_ellblock_t * eigenvectors);

void bml_diagonalize_ellblock_single_complex(
    bml_matrix_ellblock_t * A,
    void *eigenvalues,
    bml_matrix_ellblock_t * eigenvectors);

void bml_diagonalize_ellblock_double_complex(
    bml_matrix_ellblock_t * A,
    void *eigenvalues,
    bml_matrix_ellblock_t * eigenvectors);

#endif
