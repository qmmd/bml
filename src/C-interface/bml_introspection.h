/** \file */

#ifndef __BML_INTROSPECTION_H
#define __BML_INTROSPECTION_H

#include "bml_types.h"

bml_matrix_type_t bml_get_type(
    bml_matrix_t * A);

bml_matrix_precision_t bml_get_precision(
    bml_matrix_t * A);

int bml_get_N(
    bml_matrix_t * A);

int bml_get_M(
    bml_matrix_t * A);

int bml_get_NB(
    bml_matrix_t * A);

int bml_get_row_bandwidth(
    bml_matrix_t * A,
    int i);

int bml_get_bandwidth(
    bml_matrix_t * A);

double bml_get_sparsity(
    bml_matrix_t * A,
    double threshold);

bml_distribution_mode_t bml_get_distribution_mode(
    bml_matrix_t * A);

#endif
