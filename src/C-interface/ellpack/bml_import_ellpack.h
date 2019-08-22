#ifndef __BML_IMPORT_ELLPACK_H
#define __BML_IMPORT_ELLPACK_H

#include "bml_types_ellpack.h"

bml_matrix_ellpack_t *bml_import_from_dense_ellpack(
    bml_matrix_precision_t matrix_precision,
    bml_dense_order_t order,
    int N,
    void *A,
    double threshold,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_ellpack_t *bml_import_from_dense_ellpack_single_real(
    bml_dense_order_t order,
    int N,
    void *A,
    double threshold,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_ellpack_t *bml_import_from_dense_ellpack_double_real(
    bml_dense_order_t order,
    int N,
    void *A,
    double threshold,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_ellpack_t *bml_import_from_dense_ellpack_single_complex(
    bml_dense_order_t order,
    int N,
    void *A,
    double threshold,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_ellpack_t *bml_import_from_dense_ellpack_double_complex(
    bml_dense_order_t order,
    int N,
    void *A,
    double threshold,
    int M,
    bml_distribution_mode_t distrib_mode);

#endif
