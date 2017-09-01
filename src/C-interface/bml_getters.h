/** \file */

#ifndef __BML_GETTERS_H
#define __BML_GETTERS_H

#include "bml_types.h"

void *bml_get(
    const bml_matrix_t * A,
    const int i,
    const int j);

void *bml_get_row(
    bml_matrix_t * A,
    const int i);

void *bml_get_diagonal(
    bml_matrix_t * A);

#endif
