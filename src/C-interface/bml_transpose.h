/** \file */

#ifndef __BML_TRANSPOSE_H
#define __BML_TRANSPOSE_H

#include "bml_types.h"

bml_matrix_t *bml_transpose_new(
    bml_matrix_t * A);

void bml_transpose(
    bml_matrix_t * A);

void bml_complex_conjugate(bml_matrix_t * A);

#endif
