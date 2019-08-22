#ifndef __BML_INTROSPECTION_DENSE_H
#define __BML_INTROSPECTION_DENSE_H

#include "bml_types_dense.h"

bml_matrix_precision_t bml_get_precision_dense(
    bml_matrix_dense_t * A);

int bml_get_N_dense(
    bml_matrix_dense_t * A);

int bml_get_M_dense(
    bml_matrix_dense_t * A);

int bml_get_row_bandwidth_dense(
    bml_matrix_dense_t * A,
    int i);

bml_distribution_mode_t bml_get_distribution_mode_dense(
    bml_matrix_dense_t * A);

int bml_get_row_bandwidth_dense_single_real(
    bml_matrix_dense_t * A,
    int i);

int bml_get_row_bandwidth_dense_double_real(
    bml_matrix_dense_t * A,
    int i);

int bml_get_row_bandwidth_dense_single_complex(
    bml_matrix_dense_t * A,
    int i);

int bml_get_row_bandwidth_dense_double_complex(
    bml_matrix_dense_t * A,
    int i);

int bml_get_bandwidth_dense(
    bml_matrix_dense_t * A);

int bml_get_bandwidth_dense_single_real(
    bml_matrix_dense_t * A);

int bml_get_bandwidth_dense_double_real(
    bml_matrix_dense_t * A);

int bml_get_bandwidth_dense_single_complex(
    bml_matrix_dense_t * A);

int bml_get_bandwidth_dense_double_complex(
    bml_matrix_dense_t * A);

// Get the sparsity of a bml matrix
double bml_get_sparsity_dense(
    bml_matrix_dense_t * A,
    double threshold);

double bml_get_sparsity_dense_single_real(
    bml_matrix_dense_t * A,
    double threshold);

double bml_get_sparsity_dense_double_real(
    bml_matrix_dense_t * A,
    double threshold);

double bml_get_sparsity_dense_single_complex(
    bml_matrix_dense_t * A,
    double threshold);

double bml_get_sparsity_dense_double_complex(
    bml_matrix_dense_t * A,
    double threshold);

#endif
