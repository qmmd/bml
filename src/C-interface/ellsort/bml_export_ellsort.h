#ifndef __BML_EXPORT_ELLSORT_H
#define __BML_EXPORT_ELLSORT_H

#include "bml_types_ellsort.h"

void *bml_export_to_dense_ellsort(
    bml_matrix_ellsort_t * A,
    bml_dense_order_t order);

void *bml_export_to_dense_ellsort_single_real(
    bml_matrix_ellsort_t * A,
    bml_dense_order_t order);

void *bml_export_to_dense_ellsort_double_real(
    bml_matrix_ellsort_t * A,
    bml_dense_order_t order);

void *bml_export_to_dense_ellsort_single_complex(
    bml_matrix_ellsort_t * A,
    bml_dense_order_t order);

void *bml_export_to_dense_ellsort_double_complex(
    bml_matrix_ellsort_t * A,
    bml_dense_order_t order);

#endif
