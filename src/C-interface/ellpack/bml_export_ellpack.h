#ifndef __BML_EXPORT_ELLPACK_H
#define __BML_EXPORT_ELLPACK_H

#include "bml_types_ellpack.h"

void *bml_export_to_dense_ellpack(
    bml_matrix_ellpack_t * A,
    bml_dense_order_t order);

void *bml_export_to_dense_ellpack_single_real(
    bml_matrix_ellpack_t * A,
    bml_dense_order_t order);

void *bml_export_to_dense_ellpack_double_real(
    bml_matrix_ellpack_t * A,
    bml_dense_order_t order);

void *bml_export_to_dense_ellpack_single_complex(
    bml_matrix_ellpack_t * A,
    bml_dense_order_t order);

void *bml_export_to_dense_ellpack_double_complex(
    bml_matrix_ellpack_t * A,
    bml_dense_order_t order);

#endif
