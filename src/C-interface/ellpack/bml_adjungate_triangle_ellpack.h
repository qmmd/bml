/** \file */

#ifndef __BML_ADJUNGATE_TRIANGLE_ELLPACK_H
#define __BML_ADJUNGATE_TRIANGLE_ELLPACK_H

#include "bml_types_ellpack.h"

void bml_adjungate_triangle_ellpack(
    bml_matrix_ellpack_t * const A,
    const char *const triangle);

void bml_adjungate_triangle_ellpack_single_real(
    bml_matrix_ellpack_t * const A,
    const char *const triangle);

void bml_adjungate_triangle_ellpack_double_real(
    bml_matrix_ellpack_t * const A,
    const char *const triangle);

void bml_adjungate_triangle_ellpack_single_complex(
    bml_matrix_ellpack_t * const A,
    const char *const triangle);

void bml_adjungate_triangle_ellpack_double_complex(
    bml_matrix_ellpack_t * const A,
    const char *const triangle);

#endif
