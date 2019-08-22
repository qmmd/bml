/** \file */

#ifndef __BML_SCALE_H
#define __BML_SCALE_H

#include "bml_types.h"

bml_matrix_t *bml_scale_new(
    void *scale_factor,
    bml_matrix_t * A);

void bml_scale(
    void *scale_factor,
    bml_matrix_t * A,
    bml_matrix_t * B);

void bml_scale_inplace(
    void *scale_factor,
    bml_matrix_t * A);

#endif
