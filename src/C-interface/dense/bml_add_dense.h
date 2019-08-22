#ifndef __BML_ADD_DENSE_H
#define __BML_ADD_DENSE_H

#include "bml_types_dense.h"

void bml_add_dense(
    bml_matrix_dense_t * A,
    bml_matrix_dense_t * B,
    double alpha,
    double beta);

void bml_add_dense_single_real(
    bml_matrix_dense_t * A,
    bml_matrix_dense_t * B,
    double alpha,
    double beta);

void bml_add_dense_double_real(
    bml_matrix_dense_t * A,
    bml_matrix_dense_t * B,
    double alpha,
    double beta);

void bml_add_dense_single_complex(
    bml_matrix_dense_t * A,
    bml_matrix_dense_t * B,
    double alpha,
    double beta);

void bml_add_dense_double_complex(
    bml_matrix_dense_t * A,
    bml_matrix_dense_t * B,
    double alpha,
    double beta);

double bml_add_norm_dense(
    bml_matrix_dense_t * A,
    bml_matrix_dense_t * B,
    double alpha,
    double beta);

double bml_add_norm_dense_single_real(
    bml_matrix_dense_t * A,
    bml_matrix_dense_t * B,
    double alpha,
    double beta);

double bml_add_norm_dense_double_real(
    bml_matrix_dense_t * A,
    bml_matrix_dense_t * B,
    double alpha,
    double beta);

double bml_add_norm_dense_single_complex(
    bml_matrix_dense_t * A,
    bml_matrix_dense_t * B,
    double alpha,
    double beta);

double bml_add_norm_dense_double_complex(
    bml_matrix_dense_t * A,
    bml_matrix_dense_t * B,
    double alpha,
    double beta);

void bml_add_identity_dense(
    bml_matrix_dense_t * A,
    double beta);

void bml_add_identity_dense_single_real(
    bml_matrix_dense_t * A,
    double beta);

void bml_add_identity_dense_double_real(
    bml_matrix_dense_t * A,
    double beta);

void bml_add_identity_dense_single_complex(
    bml_matrix_dense_t * A,
    double beta);

void bml_add_identity_dense_double_complex(
    bml_matrix_dense_t * A,
    double beta);

void bml_scale_add_identity_dense(
    bml_matrix_dense_t * A,
    double alpha,
    double beta);

void bml_scale_add_identity_dense_single_real(
    bml_matrix_dense_t * A,
    double alpha,
    double beta);

void bml_scale_add_identity_dense_double_real(
    bml_matrix_dense_t * A,
    double alpha,
    double beta);

void bml_scale_add_identity_dense_single_complex(
    bml_matrix_dense_t * A,
    double alpha,
    double beta);

void bml_scale_add_identity_dense_double_complex(
    bml_matrix_dense_t * A,
    double alpha,
    double beta);

#endif
