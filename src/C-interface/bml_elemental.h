#ifndef __BML_ELEMENTAL
#define __BML_ELEMENTAL

#include "bml_types.h"

#include <complex.h>

float bml_get_single_real(
    bml_matrix_t * A,
    int i,
    int j);

double bml_get_double_real(
    bml_matrix_t * A,
    int i,
    int j);

float complex bml_get_single_complex(
    bml_matrix_t * A,
    int i,
    int j);

double complex bml_get_double_complex(
    bml_matrix_t * A,
    int i,
    int j);

#endif
