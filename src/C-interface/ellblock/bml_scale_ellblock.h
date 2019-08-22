#ifndef __BML_SCALE_ELLBLOCK_H
#define __BML_SCALE_ELLBLOCK_H

#include "bml_types_ellblock.h"

#include <complex.h>

bml_matrix_ellblock_t *bml_scale_ellblock_new(
    void *scale_factor,
    bml_matrix_ellblock_t * A);

bml_matrix_ellblock_t *bml_scale_ellblock_new_single_real(
    void *scale_factor,
    bml_matrix_ellblock_t * A);

bml_matrix_ellblock_t *bml_scale_ellblock_new_double_real(
    void *scale_factor,
    bml_matrix_ellblock_t * A);

bml_matrix_ellblock_t *bml_scale_ellblock_new_single_complex(
    void *scale_factor,
    bml_matrix_ellblock_t * A);

bml_matrix_ellblock_t *bml_scale_ellblock_new_double_complex(
    void *scale_factor,
    bml_matrix_ellblock_t * A);

void bml_scale_ellblock(
    void *scale_factor,
    bml_matrix_ellblock_t * A,
    bml_matrix_ellblock_t * B);

void bml_scale_ellblock_single_real(
    void *scale_factor,
    bml_matrix_ellblock_t * A,
    bml_matrix_ellblock_t * B);

void bml_scale_ellblock_double_real(
    void *scale_factor,
    bml_matrix_ellblock_t * A,
    bml_matrix_ellblock_t * B);

void bml_scale_ellblock_single_complex(
    void *scale_factor,
    bml_matrix_ellblock_t * A,
    bml_matrix_ellblock_t * B);

void bml_scale_ellblock_double_complex(
    void *scale_factor,
    bml_matrix_ellblock_t * A,
    bml_matrix_ellblock_t * B);

void bml_scale_inplace_ellblock(
    void *scale_factor,
    bml_matrix_ellblock_t * A);

void bml_scale_inplace_ellblock_single_real(
    void *scale_factor,
    bml_matrix_ellblock_t * A);

void bml_scale_inplace_ellblock_double_real(
    void *scale_factor,
    bml_matrix_ellblock_t * A);

void bml_scale_inplace_ellblock_single_complex(
    void *scale_factor,
    bml_matrix_ellblock_t * A);

void bml_scale_inplace_ellblock_double_complex(
    void *scale_factor,
    bml_matrix_ellblock_t * A);

#endif
