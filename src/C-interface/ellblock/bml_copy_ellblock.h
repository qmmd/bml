#ifndef __BML_COPY_ELLBLOCK_H
#define __BML_COPY_ELLBLOCK_H

#include "bml_types_ellblock.h"

bml_matrix_ellblock_t *bml_copy_ellblock_new(
    bml_matrix_ellblock_t * A);

bml_matrix_ellblock_t *bml_copy_ellblock_new_single_real(
    bml_matrix_ellblock_t * A);

bml_matrix_ellblock_t *bml_copy_ellblock_new_double_real(
    bml_matrix_ellblock_t * A);

bml_matrix_ellblock_t *bml_copy_ellblock_new_single_complex(
    bml_matrix_ellblock_t * A);

bml_matrix_ellblock_t *bml_copy_ellblock_new_double_complex(
    bml_matrix_ellblock_t * A);

void bml_copy_ellblock(
    bml_matrix_ellblock_t * A,
    bml_matrix_ellblock_t * B);

void bml_copy_ellblock_single_real(
    bml_matrix_ellblock_t * A,
    bml_matrix_ellblock_t * B);

void bml_copy_ellblock_double_real(
    bml_matrix_ellblock_t * A,
    bml_matrix_ellblock_t * B);

void bml_copy_ellblock_single_complex(
    bml_matrix_ellblock_t * A,
    bml_matrix_ellblock_t * B);

void bml_copy_ellblock_double_complex(
    bml_matrix_ellblock_t * A,
    bml_matrix_ellblock_t * B);

void bml_reorder_ellblock(
    bml_matrix_ellblock_t * A,
    int *perm);

void bml_reorder_ellblock_single_real(
    bml_matrix_ellblock_t * A,
    int *perm);

void bml_reorder_ellblock_double_real(
    bml_matrix_ellblock_t * A,
    int *perm);

void bml_reorder_ellblock_single_complex(
    bml_matrix_ellblock_t * A,
    int *perm);

void bml_reorder_ellblock_double_complex(
    bml_matrix_ellblock_t * A,
    int *perm);

void bml_save_domain_ellblock(
    bml_matrix_ellblock_t * A);

void bml_restore_domain_ellblock(
    bml_matrix_ellblock_t * A);

#endif
