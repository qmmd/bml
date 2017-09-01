/** \file */

#ifndef __BML_CONVERT_DENSE_H
#define __BML_CONVERT_DENSE_H

#include "bml_types_dense.h"

bml_matrix_dense_t *bml_convert_dense(
    const bml_matrix_t * A,
    const bml_matrix_precision_t matrix_precision,
    const bml_distribution_mode_t distrib_mode);

bml_matrix_dense_t *bml_convert_dense_single_real(
    const bml_matrix_t * A,
    const bml_distribution_mode_t distrib_mode);

bml_matrix_dense_t *bml_convert_dense_double_real(
    const bml_matrix_t * A,
    const bml_distribution_mode_t distrib_mode);

bml_matrix_dense_t *bml_convert_dense_single_complex(
    const bml_matrix_t * A,
    const bml_distribution_mode_t distrib_mode);

bml_matrix_dense_t *bml_convert_dense_double_complex(
    const bml_matrix_t * A,
    const bml_distribution_mode_t distrib_mode);

#endif
