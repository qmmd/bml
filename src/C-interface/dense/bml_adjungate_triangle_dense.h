#ifndef __BML_ADJUNGATE_TRIANGLE_DENSE_H
#define __BML_ADJUNGATE_TRIANGLE_DENSE_H

#include "bml_types_dense.h"

void bml_adjungate_triangle_dense(
    bml_matrix_dense_t * const A,
    const char *const triangle);

void bml_adjungate_triangle_dense_single_real(
    bml_matrix_dense_t * const A,
    const char *const triangle);

void bml_adjungate_triangle_dense_double_real(
    bml_matrix_dense_t * const A,
    const char *const triangle);

void bml_adjungate_triangle_dense_single_complex(
    bml_matrix_dense_t * const A,
    const char *const triangle);

void bml_adjungate_triangle_dense_double_complex(
    bml_matrix_dense_t * const A,
    const char *const triangle);

#endif
