/** \file */

#ifndef __BML_UTILITIES_DENSE_H
#define __BML_UTILITIES_DENSE_H

#include "bml_types_dense.h"

#include <stdio.h>
#include <stdlib.h>

void bml_print_bml_matrix_dense(
    bml_matrix_dense_t * A,
    int i_l,
    int i_u,
    int j_l,
    int j_u);

void bml_print_bml_matrix_dense_single_real(
    bml_matrix_dense_t * A,
    int i_l,
    int i_u,
    int j_l,
    int j_u);

void bml_print_bml_matrix_dense_double_real(
    bml_matrix_dense_t * A,
    int i_l,
    int i_u,
    int j_l,
    int j_u);

void bml_print_bml_matrix_dense_single_complex(
    bml_matrix_dense_t * A,
    int i_l,
    int i_u,
    int j_l,
    int j_u);

void bml_print_bml_matrix_dense_double_complex(
    bml_matrix_dense_t * A,
    int i_l,
    int i_u,
    int j_l,
    int j_u);

void bml_read_bml_matrix_dense(
    bml_matrix_dense_t * A,
    char *filename);

void bml_read_bml_matrix_dense_single_real(
    bml_matrix_dense_t * A,
    char *filename);

void bml_read_bml_matrix_dense_double_real(
    bml_matrix_dense_t * A,
    char *filename);

void bml_read_bml_matrix_dense_single_complex(
    bml_matrix_dense_t * A,
    char *filename);

void bml_read_bml_matrix_dense_double_complex(
    bml_matrix_dense_t * A,
    char *filename);

void bml_write_bml_matrix_dense(
    bml_matrix_dense_t * A,
    char *filename);

void bml_write_bml_matrix_dense_single_real(
    bml_matrix_dense_t * A,
    char *filename);

void bml_write_bml_matrix_dense_double_real(
    bml_matrix_dense_t * A,
    char *filename);

void bml_write_bml_matrix_dense_single_complex(
    bml_matrix_dense_t * A,
    char *filename);

void bml_write_bml_matrix_dense_double_complex(
    bml_matrix_dense_t * A,
    char *filename);

#endif
