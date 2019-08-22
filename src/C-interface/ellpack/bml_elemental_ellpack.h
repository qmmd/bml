#ifndef __BML_ELEMENTAL_ELLPACK_H
#define __BML_ELEMENTAL_ELLPACK_H

#include "bml_types_ellpack.h"

#include <complex.h>

float bml_get_ellpack_single_real(
    bml_matrix_ellpack_t * A,
    int i,
    int j);

double bml_get_ellpack_double_real(
    bml_matrix_ellpack_t * A,
    int i,
    int j);

float complex bml_get_ellpack_single_complex(
    bml_matrix_ellpack_t * A,
    int i,
    int j);

double complex bml_get_ellpack_double_complex(
    bml_matrix_ellpack_t * A,
    int i,
    int j);

#endif
