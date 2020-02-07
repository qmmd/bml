#ifndef __BML_IMPORT_ELLSORT_H
#define __BML_IMPORT_ELLSORT_H

#include "bml_types_ellsort.h"

bml_matrix_ellsort_t *bml_import_from_dense_ellsort(
    bml_matrix_precision_t matrix_precision,
    bml_dense_order_t order,
    int N,
    void *A,
    double threshold,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_ellsort_t
    * bml_import_from_dense_ellsort_single_real(bml_dense_order_t order,
                                                int N,
                                                void *A,
                                                double threshold,
                                                int M,
                                                bml_distribution_mode_t
                                                distrib_mode);

bml_matrix_ellsort_t
    * bml_import_from_dense_ellsort_double_real(bml_dense_order_t order,
                                                int N,
                                                void *A,
                                                double threshold,
                                                int M,
                                                bml_distribution_mode_t
                                                distrib_mode);

bml_matrix_ellsort_t
    * bml_import_from_dense_ellsort_single_complex(bml_dense_order_t order,
                                                   int N,
                                                   void *A,
                                                   double threshold,
                                                   int M,
                                                   bml_distribution_mode_t
                                                   distrib_mode);

bml_matrix_ellsort_t
    * bml_import_from_dense_ellsort_double_complex(bml_dense_order_t order,
                                                   int N,
                                                   void *A,
                                                   double threshold,
                                                   int M,
                                                   bml_distribution_mode_t
                                                   distrib_mode);

#endif
