#ifndef __BML_COPY_ELLSORT_H
#define __BML_COPY_ELLSORT_H

#include "bml_types_ellsort.h"

bml_matrix_ellsort_t *bml_copy_ellsort_new(
    bml_matrix_ellsort_t * A);

bml_matrix_ellsort_t *bml_copy_ellsort_new_single_real(
    bml_matrix_ellsort_t * A);

bml_matrix_ellsort_t *bml_copy_ellsort_new_double_real(
    bml_matrix_ellsort_t * A);

bml_matrix_ellsort_t *bml_copy_ellsort_new_single_complex(
    bml_matrix_ellsort_t * A);

bml_matrix_ellsort_t *bml_copy_ellsort_new_double_complex(
    bml_matrix_ellsort_t * A);

void bml_copy_ellsort(
    bml_matrix_ellsort_t * A,
    bml_matrix_ellsort_t * B);

void bml_copy_ellsort_single_real(
    bml_matrix_ellsort_t * A,
    bml_matrix_ellsort_t * B);

void bml_copy_ellsort_double_real(
    bml_matrix_ellsort_t * A,
    bml_matrix_ellsort_t * B);

void bml_copy_ellsort_single_complex(
    bml_matrix_ellsort_t * A,
    bml_matrix_ellsort_t * B);

void bml_copy_ellsort_double_complex(
    bml_matrix_ellsort_t * A,
    bml_matrix_ellsort_t * B);

void bml_reorder_ellsort(
    bml_matrix_ellsort_t * A,
    int *perm);

void bml_reorder_ellsort_single_real(
    bml_matrix_ellsort_t * A,
    int *perm);

void bml_reorder_ellsort_double_real(
    bml_matrix_ellsort_t * A,
    int *perm);

void bml_reorder_ellsort_single_complex(
    bml_matrix_ellsort_t * A,
    int *perm);

void bml_reorder_ellsort_double_complex(
    bml_matrix_ellsort_t * A,
    int *perm);

void bml_save_domain_ellsort(
    bml_matrix_ellsort_t * A);

void bml_restore_domain_ellsort(
    bml_matrix_ellsort_t * A);

#endif
