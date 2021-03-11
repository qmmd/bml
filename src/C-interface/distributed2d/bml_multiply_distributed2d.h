#ifndef __BML_MULTIPLY_DISTRIBUTED2D_H
#define __BML_MULTIPLY_DISTRIBUTED2D_H

#include "bml_types_distributed2d.h"

void bml_multiply_distributed2d(
    bml_matrix_distributed2d_t * A,
    bml_matrix_distributed2d_t * B,
    bml_matrix_distributed2d_t * C,
    double alpha,
    double beta,
    double threshold);

void bml_multiply_distributed2d_single_real(
    bml_matrix_distributed2d_t * A,
    bml_matrix_distributed2d_t * B,
    bml_matrix_distributed2d_t * C,
    double alpha,
    double beta,
    double threshold);

void bml_multiply_distributed2d_double_real(
    bml_matrix_distributed2d_t * A,
    bml_matrix_distributed2d_t * B,
    bml_matrix_distributed2d_t * C,
    double alpha,
    double beta,
    double threshold);

void bml_multiply_distributed2d_single_complex(
    bml_matrix_distributed2d_t * A,
    bml_matrix_distributed2d_t * B,
    bml_matrix_distributed2d_t * C,
    double alpha,
    double beta,
    double threshold);

void bml_multiply_distributed2d_double_complex(
    bml_matrix_distributed2d_t * A,
    bml_matrix_distributed2d_t * B,
    bml_matrix_distributed2d_t * C,
    double alpha,
    double beta,
    double threshold);

void *bml_multiply_x2_distributed2d(
    bml_matrix_distributed2d_t * X,
    bml_matrix_distributed2d_t * X2,
    double threshold);

void *bml_multiply_x2_distributed2d_single_real(
    bml_matrix_distributed2d_t * X,
    bml_matrix_distributed2d_t * X2,
    double threshold);

void *bml_multiply_x2_distributed2d_double_real(
    bml_matrix_distributed2d_t * X,
    bml_matrix_distributed2d_t * X2,
    double threshold);

void *bml_multiply_x2_distributed2d_single_complex(
    bml_matrix_distributed2d_t * X,
    bml_matrix_distributed2d_t * X2,
    double threshold);

void *bml_multiply_x2_distributed2d_double_complex(
    bml_matrix_distributed2d_t * X,
    bml_matrix_distributed2d_t * X2,
    double threshold);

void bml_multiply_AB_distributed2d(
    bml_matrix_distributed2d_t * A,
    bml_matrix_distributed2d_t * B,
    bml_matrix_distributed2d_t * C,
    double threshold);

void bml_multiply_AB_distributed2d_single_real(
    bml_matrix_distributed2d_t * A,
    bml_matrix_distributed2d_t * B,
    bml_matrix_distributed2d_t * C,
    double threshold);

void bml_multiply_AB_distributed2d_double_real(
    bml_matrix_distributed2d_t * A,
    bml_matrix_distributed2d_t * B,
    bml_matrix_distributed2d_t * C,
    double threshold);

void bml_multiply_AB_distributed2d_single_complex(
    bml_matrix_distributed2d_t * A,
    bml_matrix_distributed2d_t * B,
    bml_matrix_distributed2d_t * C,
    double threshold);

void bml_multiply_AB_distributed2d_double_complex(
    bml_matrix_distributed2d_t * A,
    bml_matrix_distributed2d_t * B,
    bml_matrix_distributed2d_t * C,
    double threshold);

void bml_multiply_adjust_AB_distributed2d(
    bml_matrix_distributed2d_t * A,
    bml_matrix_distributed2d_t * B,
    bml_matrix_distributed2d_t * C,
    double threshold);

void bml_multiply_adjust_AB_distributed2d_single_real(
    bml_matrix_distributed2d_t * A,
    bml_matrix_distributed2d_t * B,
    bml_matrix_distributed2d_t * C,
    double threshold);

void bml_multiply_adjust_AB_distributed2d_double_real(
    bml_matrix_distributed2d_t * A,
    bml_matrix_distributed2d_t * B,
    bml_matrix_distributed2d_t * C,
    double threshold);

void bml_multiply_adjust_AB_distributed2d_single_complex(
    bml_matrix_distributed2d_t * A,
    bml_matrix_distributed2d_t * B,
    bml_matrix_distributed2d_t * C,
    double threshold);

void bml_multiply_adjust_AB_distributed2d_double_complex(
    bml_matrix_distributed2d_t * A,
    bml_matrix_distributed2d_t * B,
    bml_matrix_distributed2d_t * C,
    double threshold);

#endif
