#ifndef __BML_MULTIPLY_ELLBLOCK_H
#define __BML_MULTIPLY_ELLBLOCK_H

#include "bml_types_ellblock.h"

void bml_multiply_ellblock(
    bml_matrix_ellblock_t * A,
    bml_matrix_ellblock_t * B,
    bml_matrix_ellblock_t * C,
    double alpha,
    double beta,
    double threshold);

void bml_multiply_ellblock_single_real(
    bml_matrix_ellblock_t * A,
    bml_matrix_ellblock_t * B,
    bml_matrix_ellblock_t * C,
    double alpha,
    double beta,
    double threshold);

void bml_multiply_ellblock_double_real(
    bml_matrix_ellblock_t * A,
    bml_matrix_ellblock_t * B,
    bml_matrix_ellblock_t * C,
    double alpha,
    double beta,
    double threshold);

void bml_multiply_ellblock_single_complex(
    bml_matrix_ellblock_t * A,
    bml_matrix_ellblock_t * B,
    bml_matrix_ellblock_t * C,
    double alpha,
    double beta,
    double threshold);

void bml_multiply_ellblock_double_complex(
    bml_matrix_ellblock_t * A,
    bml_matrix_ellblock_t * B,
    bml_matrix_ellblock_t * C,
    double alpha,
    double beta,
    double threshold);

void *bml_multiply_x2_ellblock(
    bml_matrix_ellblock_t * X,
    bml_matrix_ellblock_t * X2,
    double threshold);

void *bml_multiply_x2_ellblock_single_real(
    bml_matrix_ellblock_t * X,
    bml_matrix_ellblock_t * X2,
    double threshold);

void *bml_multiply_x2_ellblock_double_real(
    bml_matrix_ellblock_t * X,
    bml_matrix_ellblock_t * X2,
    double threshold);

void *bml_multiply_x2_ellblock_single_complex(
    bml_matrix_ellblock_t * X,
    bml_matrix_ellblock_t * X2,
    double threshold);

void *bml_multiply_x2_ellblock_double_complex(
    bml_matrix_ellblock_t * X,
    bml_matrix_ellblock_t * X2,
    double threshold);

void bml_multiply_AB_ellblock(
    bml_matrix_ellblock_t * A,
    bml_matrix_ellblock_t * B,
    bml_matrix_ellblock_t * C,
    double threshold);

void bml_multiply_AB_ellblock_single_real(
    bml_matrix_ellblock_t * A,
    bml_matrix_ellblock_t * B,
    bml_matrix_ellblock_t * C,
    double threshold);

void bml_multiply_AB_ellblock_double_real(
    bml_matrix_ellblock_t * A,
    bml_matrix_ellblock_t * B,
    bml_matrix_ellblock_t * C,
    double threshold);

void bml_multiply_AB_ellblock_single_complex(
    bml_matrix_ellblock_t * A,
    bml_matrix_ellblock_t * B,
    bml_matrix_ellblock_t * C,
    double threshold);

void bml_multiply_AB_ellblock_double_complex(
    bml_matrix_ellblock_t * A,
    bml_matrix_ellblock_t * B,
    bml_matrix_ellblock_t * C,
    double threshold);

void bml_multiply_adjust_AB_ellblock(
    bml_matrix_ellblock_t * A,
    bml_matrix_ellblock_t * B,
    bml_matrix_ellblock_t * C,
    double threshold);

void bml_multiply_adjust_AB_ellblock_single_real(
    bml_matrix_ellblock_t * A,
    bml_matrix_ellblock_t * B,
    bml_matrix_ellblock_t * C,
    double threshold);

void bml_multiply_adjust_AB_ellblock_double_real(
    bml_matrix_ellblock_t * A,
    bml_matrix_ellblock_t * B,
    bml_matrix_ellblock_t * C,
    double threshold);

void bml_multiply_adjust_AB_ellblock_single_complex(
    bml_matrix_ellblock_t * A,
    bml_matrix_ellblock_t * B,
    bml_matrix_ellblock_t * C,
    double threshold);

void bml_multiply_adjust_AB_ellblock_double_complex(
    bml_matrix_ellblock_t * A,
    bml_matrix_ellblock_t * B,
    bml_matrix_ellblock_t * C,
    double threshold);

#endif
