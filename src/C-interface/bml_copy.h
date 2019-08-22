/** \file */

#ifndef __BML_COPY_H
#define __BML_COPY_H

#include "bml_types.h"

bml_matrix_t *bml_copy_new(
    bml_matrix_t * A);

void bml_copy(
    bml_matrix_t * A,
    bml_matrix_t * B);

void bml_reorder(
    bml_matrix_t * A,
    int *perm);

void bml_copy_domain(
    bml_domain_t * A,
    bml_domain_t * B);

void bml_save_domain(
    bml_matrix_t * A);

void bml_restore_domain(
    bml_matrix_t * A);

#endif
