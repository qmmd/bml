#ifndef __BML_IMPORT_DENSE_H
#define __BML_IMPORT_DENSE_H

#include "bml_types_dense.h"

bml_matrix_dense_t *bml_import_from_dense_dense(
    bml_matrix_precision_t matrix_precision,
    bml_dense_order_t order,
    int N,
    void *A,
    double threshold,
    bml_distribution_mode_t distrib_mode);

bml_matrix_dense_t *bml_import_from_dense_dense_single_real(
    bml_dense_order_t order,
    int N,
    void *A,
    bml_distribution_mode_t distrib_mode);

bml_matrix_dense_t *bml_import_from_dense_dense_double_real(
    bml_dense_order_t order,
    int N,
    void *A,
    bml_distribution_mode_t distrib_mode);

bml_matrix_dense_t *bml_import_from_dense_dense_single_complex(
    bml_dense_order_t order,
    int N,
    void *A,
    bml_distribution_mode_t distrib_mode);

bml_matrix_dense_t *bml_import_from_dense_dense_double_complex(
    bml_dense_order_t order,
    int N,
    void *A,
    bml_distribution_mode_t distrib_mode);

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
