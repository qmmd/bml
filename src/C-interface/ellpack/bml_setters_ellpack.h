/** \file */

#ifndef __BML_SETTERS_ELLPACK_H
#define __BML_SETTERS_ELLPACK_H

#include "bml_types_ellpack.h"

#include <complex.h>

void bml_set_element_new_ellpack(
    bml_matrix_ellpack_t * A,
    int i,
    int j,
    void *value);

void bml_set_element_new_ellpack_single_real(
    bml_matrix_ellpack_t * A,
    int i,
    int j,
    void *value);

void bml_set_element_new_ellpack_double_real(
    bml_matrix_ellpack_t * A,
    int i,
    int j,
    void *value);

void bml_set_element_new_ellpack_single_complex(
    bml_matrix_ellpack_t * A,
    int i,
    int j,
    void *value);

void bml_set_element_new_ellpack_double_complex(
    bml_matrix_ellpack_t * A,
    int i,
    int j,
    void *value);

void bml_set_element_ellpack(
    bml_matrix_ellpack_t * A,
    int i,
    int j,
    void *value);

void bml_set_element_ellpack_single_real(
    bml_matrix_ellpack_t * A,
    int i,
    int j,
    void *value);

void bml_set_element_ellpack_double_real(
    bml_matrix_ellpack_t * A,
    int i,
    int j,
    void *value);

void bml_set_element_ellpack_single_complex(
    bml_matrix_ellpack_t * A,
    int i,
    int j,
    void *value);

void bml_set_element_ellpack_double_complex(
    bml_matrix_ellpack_t * A,
    int i,
    int j,
    void *value);

void bml_set_row_ellpack(
    bml_matrix_ellpack_t * A,
    int i,
    void *row,
    double threshold);

void bml_set_row_ellpack_single_real(
    bml_matrix_ellpack_t * A,
    int i,
    void *row,
    double threshold);

void bml_set_row_ellpack_double_real(
    bml_matrix_ellpack_t * A,
    int i,
    void *row,
    double threshold);

void bml_set_row_ellpack_single_complex(
    bml_matrix_ellpack_t * A,
    int i,
    void *row,
    double threshold);

void bml_set_row_ellpack_double_complex(
    bml_matrix_ellpack_t * A,
    int i,
    void *row,
    double threshold);

void bml_set_diagonal_ellpack(
    bml_matrix_ellpack_t * A,
    void *diagonal,
    double threshold);

void bml_set_diagonal_ellpack_single_real(
    bml_matrix_ellpack_t * A,
    void *diagonal,
    double threshold);

void bml_set_diagonal_ellpack_double_real(
    bml_matrix_ellpack_t * A,
    void *diagonal,
    double threshold);

void bml_set_diagonal_ellpack_single_complex(
    bml_matrix_ellpack_t * A,
    void *diagonal,
    double threshold);

void bml_set_diagonal_ellpack_double_complex(
    bml_matrix_ellpack_t * A,
    void *diagonal,
    double threshold);

#endif
