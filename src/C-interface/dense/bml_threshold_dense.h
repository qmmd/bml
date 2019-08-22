#ifndef __BML_THRESHOLD_DENSE_H
#define __BML_THRESHOLD_DENSE_H

#include "bml_types_dense.h"

bml_matrix_dense_t *bml_threshold_new_dense(
    bml_matrix_dense_t * A,
    double threshold);

bml_matrix_dense_t *bml_threshold_new_dense_single_real(
    bml_matrix_dense_t * A,
    double threshold);

bml_matrix_dense_t *bml_threshold_new_dense_double_real(
    bml_matrix_dense_t * A,
    double threshold);

bml_matrix_dense_t *bml_threshold_new_dense_single_complex(
    bml_matrix_dense_t * A,
    double threshold);

bml_matrix_dense_t *bml_threshold_new_dense_double_complex(
    bml_matrix_dense_t * A,
    double threshold);

void bml_threshold_dense(
    bml_matrix_dense_t * A,
    double threshold);

void bml_threshold_dense_single_real(
    bml_matrix_dense_t * A,
    double threshold);

void bml_threshold_dense_double_real(
    bml_matrix_dense_t * A,
    double threshold);

void bml_threshold_dense_single_complex(
    bml_matrix_dense_t * A,
    double threshold);

void bml_threshold_dense_double_complex(
    bml_matrix_dense_t * A,
    double threshold);

#endif
