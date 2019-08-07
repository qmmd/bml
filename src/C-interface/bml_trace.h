/** \file */

#ifndef __BML_TRACE_H
#define __BML_TRACE_H

#include "bml_types.h"

// Calculate trace of A
double bml_trace(
    bml_matrix_t * A);

// Calculate trace of matrix mult
double bml_trace_mult(
    bml_matrix_t * A,
    bml_matrix_t * B);

#endif
