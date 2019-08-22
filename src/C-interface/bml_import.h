/** \file */

#ifndef __BML_IMPORT_H
#define __BML_IMPORT_H

#include "bml_types.h"

bml_matrix_t *bml_import_from_dense(
    bml_matrix_type_t matrix_type,
    bml_matrix_precision_t matrix_precision,
    bml_dense_order_t order,
    int N,
    int M,
    void *A,
    double threshold,
    bml_distribution_mode_t distrib_mode);

#endif
