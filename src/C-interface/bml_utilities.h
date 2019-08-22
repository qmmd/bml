/** \file */

#ifndef __BML_UTILITIES_H
#define __BML_UTILITIES_H

#include "bml_types.h"

void bml_print_dense_matrix(
    int N,
    bml_matrix_precision_t matrix_precision,
    bml_dense_order_t order,
    void *A,
    int i_l,
    int i_u,
    int j_l,
    int j_u);

void bml_print_dense_vector(
    int N,
    bml_matrix_precision_t matrix_precision,
    void *v,
    int i_l,
    int i_u);

void bml_print_bml_vector(
    bml_vector_t * v,
    int i_l,
    int i_u);

void bml_print_bml_matrix(
    bml_matrix_t * A,
    int i_l,
    int i_u,
    int j_l,
    int j_u);

void bml_read_bml_matrix(
    bml_matrix_t * A,
    char *filename);

void bml_write_bml_matrix(
    bml_matrix_t * A,
    char *filename);

#endif
