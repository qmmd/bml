/** \file */

#ifndef __BML_CONVERT_ELLPACK_H
#define __BML_CONVERT_ELLPACK_H

#include "bml_types_ellpack.h"

bml_matrix_ellpack_t *bml_convert_ellpack(
    const bml_matrix_t * A,
    const bml_matrix_precision_t matrix_precision,
    const int M,
    const bml_distribution_mode_t distrib_mode);

bml_matrix_ellpack_t *bml_convert_ellpack_single_real(
    const bml_matrix_t * A,
    const bml_matrix_precision_t matrix_precision,
    const int M,
    const bml_distribution_mode_t distrib_mode);

bml_matrix_ellpack_t *bml_convert_ellpack_double_real(
    const bml_matrix_t * A,
    const bml_matrix_precision_t matrix_precision,
    const int M,
    const bml_distribution_mode_t distrib_mode);

bml_matrix_ellpack_t *bml_convert_ellpack_single_complex(
    const bml_matrix_t * A,
    const bml_matrix_precision_t matrix_precision,
    const int M,
    const bml_distribution_mode_t distrib_mode);

bml_matrix_ellpack_t *bml_convert_ellpack_double_complex(
    const bml_matrix_t * A,
    const bml_matrix_precision_t matrix_precision,
    const int M,
    const bml_distribution_mode_t distrib_mode);

#endif
