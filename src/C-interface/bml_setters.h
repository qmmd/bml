/** \file */

#ifndef __BML_SETTERS_H
#define __BML_SETTERS_H

#include "bml_types.h"

void bml_set_element_new(
    bml_matrix_t * A,
    int i,
    int j,
    void *value);

void bml_set_element(
    bml_matrix_t * A,
    int i,
    int j,
    void *value);

void bml_set_row(
    bml_matrix_t * A,
    int i,
    void *row,
    double threshold);

void bml_set_diagonal(
    bml_matrix_t * A,
    void *diagonal,
    double threshold);

#endif
