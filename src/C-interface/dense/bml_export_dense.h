#ifndef __BML_EXPORT_DENSE_H
#define __BML_EXPORT_DENSE_H

#include "bml_types_dense.h"

void *bml_export_to_dense_dense(
    bml_matrix_dense_t * A,
    bml_dense_order_t order);

void *bml_export_to_dense_dense_single_real(
    bml_matrix_dense_t * A,
    bml_dense_order_t order);

void *bml_export_to_dense_dense_double_real(
    bml_matrix_dense_t * A,
    bml_dense_order_t order);

void *bml_export_to_dense_dense_single_complex(
    bml_matrix_dense_t * A,
    bml_dense_order_t order);

void *bml_export_to_dense_dense_double_complex(
    bml_matrix_dense_t * A,
    bml_dense_order_t order);

#endif
