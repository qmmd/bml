/** \file */

#ifndef __BML_GETTERS_H
#define __BML_GETTERS_H

#include "bml_types.h"

void *bml_get(
    bml_matrix_t * A,
    int i,
    int j);

void *bml_get_row(
    bml_matrix_t * A,
    int i);

void *bml_get_diagonal(
    bml_matrix_t * A);

#endif
