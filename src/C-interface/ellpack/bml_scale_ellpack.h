#ifndef __BML_SCALE_ELLPACK_H
#define __BML_SCALE_ELLPACK_H

#include "bml_types_ellpack.h"

#include <complex.h>

bml_matrix_ellpack_t *bml_scale_ellpack_new(
    void *scale_factor,
    bml_matrix_ellpack_t * A);

bml_matrix_ellpack_t *bml_scale_ellpack_new_single_real(
    void *scale_factor,
    bml_matrix_ellpack_t * A);

bml_matrix_ellpack_t *bml_scale_ellpack_new_double_real(
    void *scale_factor,
    bml_matrix_ellpack_t * A);

bml_matrix_ellpack_t *bml_scale_ellpack_new_single_complex(
    void *scale_factor,
    bml_matrix_ellpack_t * A);

bml_matrix_ellpack_t *bml_scale_ellpack_new_double_complex(
    void *scale_factor,
    bml_matrix_ellpack_t * A);

void bml_scale_ellpack(
    void *scale_factor,
    bml_matrix_ellpack_t * A,
    bml_matrix_ellpack_t * B);

void bml_scale_ellpack_single_real(
    void *scale_factor,
    bml_matrix_ellpack_t * A,
    bml_matrix_ellpack_t * B);

void bml_scale_ellpack_double_real(
    void *scale_factor,
    bml_matrix_ellpack_t * A,
    bml_matrix_ellpack_t * B);

void bml_scale_ellpack_single_complex(
    void *scale_factor,
    bml_matrix_ellpack_t * A,
    bml_matrix_ellpack_t * B);

void bml_scale_ellpack_double_complex(
    void *scale_factor,
    bml_matrix_ellpack_t * A,
    bml_matrix_ellpack_t * B);

void bml_scale_inplace_ellpack(
    void *scale_factor,
    bml_matrix_ellpack_t * A);

void bml_scale_inplace_ellpack_single_real(
    void *scale_factor,
    bml_matrix_ellpack_t * A);

void bml_scale_inplace_ellpack_double_real(
    void *scale_factor,
    bml_matrix_ellpack_t * A);

void bml_scale_inplace_ellpack_single_complex(
    void *scale_factor,
    bml_matrix_ellpack_t * A);

void bml_scale_inplace_ellpack_double_complex(
    void *scale_factor,
    bml_matrix_ellpack_t * A);

#endif
