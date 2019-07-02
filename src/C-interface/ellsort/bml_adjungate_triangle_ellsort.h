/** \file */

#ifndef __BML_ADJUNGATE_TRIANGLE_ELLSORT_H
#define __BML_ADJUNGATE_TRIANGLE_ELLSORT_H

#include "bml_types_ellsort.h"

void bml_adjungate_triangle_ellsort(
    bml_matrix_ellsort_t * const A,
    const char *const triangle);

void bml_adjungate_triangle_ellsort_single_real(
    bml_matrix_ellsort_t * const A,
    const char *const triangle);

void bml_adjungate_triangle_ellsort_double_real(
    bml_matrix_ellsort_t * const A,
    const char *const triangle);

void bml_adjungate_triangle_ellsort_single_complex(
    bml_matrix_ellsort_t * const A,
    const char *const triangle);

void bml_adjungate_triangle_ellsort_double_complex(
    bml_matrix_ellsort_t * const A,
    const char *const triangle);

#endif
