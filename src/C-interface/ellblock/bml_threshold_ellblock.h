#ifndef __BML_THRESHOLD_ELLBLOCK_H
#define __BML_THRESHOLD_ELLBLOCK_H

#include "bml_types_ellblock.h"

bml_matrix_ellblock_t *bml_threshold_new_ellblock(
    bml_matrix_ellblock_t * A,
    double threshold);

bml_matrix_ellblock_t *bml_threshold_new_ellblock_single_real(
    bml_matrix_ellblock_t * A,
    double threshold);

bml_matrix_ellblock_t *bml_threshold_new_ellblock_double_real(
    bml_matrix_ellblock_t * A,
    double threshold);

bml_matrix_ellblock_t *bml_threshold_new_ellblock_single_complex(
    bml_matrix_ellblock_t * A,
    double threshold);

bml_matrix_ellblock_t *bml_threshold_new_ellblock_double_complex(
    bml_matrix_ellblock_t * A,
    double threshold);

void bml_threshold_ellblock(
    bml_matrix_ellblock_t * A,
    double threshold);

void bml_threshold_ellblock_single_real(
    bml_matrix_ellblock_t * A,
    double threshold);

void bml_threshold_ellblock_double_real(
    bml_matrix_ellblock_t * A,
    double threshold);

void bml_threshold_ellblock_single_complex(
    bml_matrix_ellblock_t * A,
    double threshold);

void bml_threshold_ellblock_double_complex(
    bml_matrix_ellblock_t * A,
    double threshold);

#endif
