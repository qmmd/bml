#ifndef __BML_TRACE_DENSE_H
#define __BML_TRACE_DENSE_H

#include "bml_types_dense.h"

double bml_trace_dense(
    bml_matrix_dense_t * A);

double bml_trace_dense_single_real(
    bml_matrix_dense_t * A);

double bml_trace_dense_double_real(
    bml_matrix_dense_t * A);

double bml_trace_dense_single_complex(
    bml_matrix_dense_t * A);

double bml_trace_dense_double_complex(
    bml_matrix_dense_t * A);

double bml_trace_mult_dense(
    bml_matrix_dense_t * A,
    bml_matrix_dense_t * B);

double bml_trace_mult_dense_single_real(
    bml_matrix_dense_t * A,
    bml_matrix_dense_t * B);

double bml_trace_mult_dense_double_real(
    bml_matrix_dense_t * A,
    bml_matrix_dense_t * B);

double bml_trace_mult_dense_single_complex(
    bml_matrix_dense_t * A,
    bml_matrix_dense_t * B);

double bml_trace_mult_dense_double_complex(
    bml_matrix_dense_t * A,
    bml_matrix_dense_t * B);

#endif
