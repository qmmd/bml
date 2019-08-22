/** \file */

#ifndef __BML_GETTERS_ELLPACK_H
#define __BML_GETTERS_ELLPACK_H

#include "bml_types_ellpack.h"

#include <complex.h>

void *bml_get_ellpack(
    bml_matrix_ellpack_t * A,
    int i,
    int j);

void *bml_get_ellpack_single_real(
    bml_matrix_ellpack_t * A,
    int i,
    int j);

void *bml_get_ellpack_double_real(
    bml_matrix_ellpack_t * A,
    int i,
    int j);

void *bml_get_ellpack_single_complex(
    bml_matrix_ellpack_t * A,
    int i,
    int j);

void *bml_get_ellpack_double_complex(
    bml_matrix_ellpack_t * A,
    int i,
    int j);

void *bml_get_row_ellpack(
    bml_matrix_ellpack_t * A,
    int i);

void *bml_get_row_ellpack_single_real(
    bml_matrix_ellpack_t * A,
    int i);

void *bml_get_row_ellpack_double_real(
    bml_matrix_ellpack_t * A,
    int i);

void *bml_get_row_ellpack_single_complex(
    bml_matrix_ellpack_t * A,
    int i);

void *bml_get_row_ellpack_double_complex(
    bml_matrix_ellpack_t * A,
    int i);

void *bml_get_diagonal_ellpack(
    bml_matrix_ellpack_t * A);

void *bml_get_diagonal_ellpack_single_real(
    bml_matrix_ellpack_t * A);

void *bml_get_diagonal_ellpack_double_real(
    bml_matrix_ellpack_t * A);

void *bml_get_diagonal_ellpack_single_complex(
    bml_matrix_ellpack_t * A);

void *bml_get_diagonal_ellpack_double_complex(
    bml_matrix_ellpack_t * A);

#endif
