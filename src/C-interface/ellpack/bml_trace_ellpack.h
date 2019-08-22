#ifndef __BML_TRACE_ELLPACK_H
#define __BML_TRACE_ELLPACK_H

#include "bml_types_ellpack.h"

double bml_trace_ellpack(
    bml_matrix_ellpack_t * A);

double bml_trace_ellpack_single_real(
    bml_matrix_ellpack_t * A);

double bml_trace_ellpack_double_real(
    bml_matrix_ellpack_t * A);

double bml_trace_ellpack_single_complex(
    bml_matrix_ellpack_t * A);

double bml_trace_ellpack_double_complex(
    bml_matrix_ellpack_t * A);

double bml_trace_mult_ellpack(
    bml_matrix_ellpack_t * A,
    bml_matrix_ellpack_t * B);

double bml_trace_mult_ellpack_single_real(
    bml_matrix_ellpack_t * A,
    bml_matrix_ellpack_t * B);

double bml_trace_mult_ellpack_double_real(
    bml_matrix_ellpack_t * A,
    bml_matrix_ellpack_t * B);

double bml_trace_mult_ellpack_single_complex(
    bml_matrix_ellpack_t * A,
    bml_matrix_ellpack_t * B);

double bml_trace_mult_ellpack_double_complex(
    bml_matrix_ellpack_t * A,
    bml_matrix_ellpack_t * B);

#endif
