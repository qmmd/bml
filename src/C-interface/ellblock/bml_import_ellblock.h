#ifndef __BML_IMPORT_ELLBLOCK_H
#define __BML_IMPORT_ELLBLOCK_H

#include "bml_types_ellblock.h"

bml_matrix_ellblock_t *bml_import_from_dense_ellblock(
    bml_matrix_precision_t matrix_precision,
    bml_dense_order_t order,
    int N,
    void *A,
    double threshold,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_ellblock_t
    * bml_import_from_dense_ellblock_single_real(bml_dense_order_t order,
                                                 int N,
                                                 void *A,
                                                 double threshold,
                                                 int M,
                                                 bml_distribution_mode_t
                                                 distrib_mode);

bml_matrix_ellblock_t
    * bml_import_from_dense_ellblock_double_real(bml_dense_order_t order,
                                                 int N,
                                                 void *A,
                                                 double threshold,
                                                 int M,
                                                 bml_distribution_mode_t
                                                 distrib_mode);

bml_matrix_ellblock_t
    * bml_import_from_dense_ellblock_single_complex(bml_dense_order_t order,
                                                    int N,
                                                    void *A,
                                                    double threshold,
                                                    int M,
                                                    bml_distribution_mode_t
                                                    distrib_mode);

bml_matrix_ellblock_t
    * bml_import_from_dense_ellblock_double_complex(bml_dense_order_t order,
                                                    int N,
                                                    void *A,
                                                    double threshold,
                                                    int M,
                                                    bml_distribution_mode_t
                                                    distrib_mode);

#endif
