#ifndef __BML_ALLOCATE_DISTRIBUTED2D_H
#define __BML_ALLOCATE_DISTRIBUTED2D_H

#include "bml_types_distributed2d.h"

void bml_setcomm_distributed2d(
    MPI_Comm comm);

void bml_setup_distributed2d(
    const int N,
    bml_matrix_distributed2d_t * A);

void bml_deallocate_distributed2d(
    bml_matrix_distributed2d_t * A);

void bml_clear_distributed2d(
    bml_matrix_distributed2d_t * A);

void bml_clear_distributed2d_single_real(
    bml_matrix_distributed2d_t * A);

void bml_clear_distributed2d_double_real(
    bml_matrix_distributed2d_t * A);

void bml_clear_distributed2d_single_complex(
    bml_matrix_distributed2d_t * A);

void bml_clear_distributed2d_double_complex(
    bml_matrix_distributed2d_t * A);

bml_matrix_distributed2d_t *bml_zero_matrix_distributed2d(
    bml_matrix_type_t matrix_type,
    bml_matrix_precision_t matrix_precision,
    int N,
    int M);

bml_matrix_distributed2d_t *bml_random_matrix_distributed2d(
    bml_matrix_type_t matrix_type,
    bml_matrix_precision_t matrix_precision,
    int N,
    int M);

bml_matrix_distributed2d_t *bml_identity_matrix_distributed2d(
    bml_matrix_type_t matrix_type,
    bml_matrix_precision_t matrix_precision,
    int N,
    int M);

#endif
