#ifndef __BML_ELEMENTAL_DENSE_H
#define __BML_ELEMENTAL_DENSE_H

#include "bml_types_dense.h"

#include <complex.h>

float bml_get_dense_single_real(
    bml_matrix_dense_t * A,
    int i,
    int j);

double bml_get_dense_double_real(
    bml_matrix_dense_t * A,
    int i,
    int j);

float complex bml_get_dense_single_complex(
    bml_matrix_dense_t * A,
    int i,
    int j);

double complex bml_get_dense_double_complex(
    bml_matrix_dense_t * A,
    int i,
    int j);

#endif
