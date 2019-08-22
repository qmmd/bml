#ifndef __BML_COPY_ELLPACK_H
#define __BML_COPY_ELLPACK_H

#include "bml_types_ellpack.h"

bml_matrix_ellpack_t *bml_copy_ellpack_new(
    bml_matrix_ellpack_t * A);

bml_matrix_ellpack_t *bml_copy_ellpack_new_single_real(
    bml_matrix_ellpack_t * A);

bml_matrix_ellpack_t *bml_copy_ellpack_new_double_real(
    bml_matrix_ellpack_t * A);

bml_matrix_ellpack_t *bml_copy_ellpack_new_single_complex(
    bml_matrix_ellpack_t * A);

bml_matrix_ellpack_t *bml_copy_ellpack_new_double_complex(
    bml_matrix_ellpack_t * A);

void bml_copy_ellpack(
    bml_matrix_ellpack_t * A,
    bml_matrix_ellpack_t * B);

void bml_copy_ellpack_single_real(
    bml_matrix_ellpack_t * A,
    bml_matrix_ellpack_t * B);

void bml_copy_ellpack_double_real(
    bml_matrix_ellpack_t * A,
    bml_matrix_ellpack_t * B);

void bml_copy_ellpack_single_complex(
    bml_matrix_ellpack_t * A,
    bml_matrix_ellpack_t * B);

void bml_copy_ellpack_double_complex(
    bml_matrix_ellpack_t * A,
    bml_matrix_ellpack_t * B);

void bml_reorder_ellpack(
    bml_matrix_ellpack_t * A,
    int *perm);

void bml_reorder_ellpack_single_real(
    bml_matrix_ellpack_t * A,
    int *perm);

void bml_reorder_ellpack_double_real(
    bml_matrix_ellpack_t * A,
    int *perm);

void bml_reorder_ellpack_single_complex(
    bml_matrix_ellpack_t * A,
    int *perm);

void bml_reorder_ellpack_double_complex(
    bml_matrix_ellpack_t * A,
    int *perm);

void bml_save_domain_ellpack(
    bml_matrix_ellpack_t * A);

void bml_restore_domain_ellpack(
    bml_matrix_ellpack_t * A);

#endif
