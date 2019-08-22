/** \file */

#ifndef __BML_CONVERT_DENSE_H
#define __BML_CONVERT_DENSE_H

#include "bml_types_dense.h"

bml_matrix_dense_t *bml_convert_dense(
    bml_matrix_t * A,
    bml_matrix_precision_t matrix_precision,
    bml_distribution_mode_t distrib_mode);

bml_matrix_dense_t *bml_convert_dense_single_real(
    bml_matrix_t * A,
    bml_matrix_precision_t matrix_precision,
    bml_distribution_mode_t distrib_mode);

bml_matrix_dense_t *bml_convert_dense_double_real(
    bml_matrix_t * A,
    bml_matrix_precision_t matrix_precision,
    bml_distribution_mode_t distrib_mode);

bml_matrix_dense_t *bml_convert_dense_single_complex(
    bml_matrix_t * A,
    bml_matrix_precision_t matrix_precision,
    bml_distribution_mode_t distrib_mode);

bml_matrix_dense_t *bml_convert_dense_double_complex(
    bml_matrix_t * A,
    bml_matrix_precision_t matrix_precision,
    bml_distribution_mode_t distrib_mode);

#endif
