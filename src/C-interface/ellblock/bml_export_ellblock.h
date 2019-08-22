#ifndef __BML_EXPORT_ELLBLOCK_H
#define __BML_EXPORT_ELLBLOCK_H

#include "bml_types_ellblock.h"

void *bml_export_to_dense_ellblock(
    bml_matrix_ellblock_t * A,
    bml_dense_order_t order);

void *bml_export_to_dense_ellblock_single_real(
    bml_matrix_ellblock_t * A,
    bml_dense_order_t order);

void *bml_export_to_dense_ellblock_double_real(
    bml_matrix_ellblock_t * A,
    bml_dense_order_t order);

void *bml_export_to_dense_ellblock_single_complex(
    bml_matrix_ellblock_t * A,
    bml_dense_order_t order);

void *bml_export_to_dense_ellblock_double_complex(
    bml_matrix_ellblock_t * A,
    bml_dense_order_t order);

#endif
