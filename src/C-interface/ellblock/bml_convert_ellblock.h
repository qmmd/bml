/** \file */

#ifndef __BML_CONVERT_ELLBLOCK_H
#define __BML_CONVERT_ELLBLOCK_H

#include "bml_types_ellblock.h"

bml_matrix_ellblock_t *bml_convert_ellblock(
    bml_matrix_t * A,
    bml_matrix_precision_t matrix_precision,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_ellblock_t *bml_convert_ellblock_single_real(
    bml_matrix_t * A,
    bml_matrix_precision_t matrix_precision,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_ellblock_t *bml_convert_ellblock_double_real(
    bml_matrix_t * A,
    bml_matrix_precision_t matrix_precision,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_ellblock_t *bml_convert_ellblock_single_complex(
    bml_matrix_t * A,
    bml_matrix_precision_t matrix_precision,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_ellblock_t *bml_convert_ellblock_double_complex(
    bml_matrix_t * A,
    bml_matrix_precision_t matrix_precision,
    int M,
    bml_distribution_mode_t distrib_mode);

#endif
