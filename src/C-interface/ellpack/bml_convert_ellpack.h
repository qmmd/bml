/** \file */

#ifndef __BML_CONVERT_ELLPACK_H
#define __BML_CONVERT_ELLPACK_H

#include "bml_types_ellpack.h"

bml_matrix_ellpack_t *bml_convert_ellpack(
    bml_matrix_t * A,
    bml_matrix_precision_t matrix_precision,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_ellpack_t *bml_convert_ellpack_single_real(
    bml_matrix_t * A,
    bml_matrix_precision_t matrix_precision,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_ellpack_t *bml_convert_ellpack_double_real(
    bml_matrix_t * A,
    bml_matrix_precision_t matrix_precision,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_ellpack_t *bml_convert_ellpack_single_complex(
    bml_matrix_t * A,
    bml_matrix_precision_t matrix_precision,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_ellpack_t *bml_convert_ellpack_double_complex(
    bml_matrix_t * A,
    bml_matrix_precision_t matrix_precision,
    int M,
    bml_distribution_mode_t distrib_mode);

#endif
