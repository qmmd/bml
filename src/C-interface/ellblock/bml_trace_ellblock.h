#ifndef __BML_TRACE_ELLBLOCK_H
#define __BML_TRACE_ELLNLOCK_H

#include "bml_types_ellblock.h"

double bml_trace_ellblock(
    bml_matrix_ellblock_t * A);

double bml_trace_ellblock_single_real(
    bml_matrix_ellblock_t * A);

double bml_trace_ellblock_double_real(
    bml_matrix_ellblock_t * A);

double bml_trace_ellblock_single_complex(
    bml_matrix_ellblock_t * A);

double bml_trace_ellblock_double_complex(
    bml_matrix_ellblock_t * A);

double bml_trace_mult_ellblock(
    bml_matrix_ellblock_t * A,
    bml_matrix_ellblock_t * B);

double bml_trace_mult_ellblock_single_real(
    bml_matrix_ellblock_t * A,
    bml_matrix_ellblock_t * B);

double bml_trace_mult_ellblock_double_real(
    bml_matrix_ellblock_t * A,
    bml_matrix_ellblock_t * B);

double bml_trace_mult_ellblock_single_complex(
    bml_matrix_ellblock_t * A,
    bml_matrix_ellblock_t * B);

double bml_trace_mult_ellblock_double_complex(
    bml_matrix_ellblock_t * A,
    bml_matrix_ellblock_t * B);

#endif
