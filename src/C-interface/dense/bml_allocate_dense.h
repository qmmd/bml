#ifndef __BML_ALLOCATE_DENSE_H
#define __BML_ALLOCATE_DENSE_H

#include "bml_types_dense.h"

void bml_deallocate_dense(
    bml_matrix_dense_t * A);

void bml_clear_dense(
    bml_matrix_dense_t * A);

void bml_clear_dense_single_real(
    bml_matrix_dense_t * A);

void bml_clear_dense_double_real(
    bml_matrix_dense_t * A);

void bml_clear_dense_single_complex(
    bml_matrix_dense_t * A);

void bml_clear_dense_double_complex(
    bml_matrix_dense_t * A);

bml_matrix_dense_t *bml_zero_matrix_dense(
    bml_matrix_precision_t matrix_precision,
    bml_matrix_dimension_t matrix_dimension,
    bml_distribution_mode_t distrib_mode);

bml_matrix_dense_t *bml_zero_matrix_dense_single_real(
    bml_matrix_dimension_t matrix_dimension,
    bml_distribution_mode_t distrib_mode);

bml_matrix_dense_t *bml_zero_matrix_dense_double_real(
    bml_matrix_dimension_t matrix_dimension,
    bml_distribution_mode_t distrib_mode);

bml_matrix_dense_t *bml_zero_matrix_dense_single_complex(
    bml_matrix_dimension_t matrix_dimension,
    bml_distribution_mode_t distrib_mode);

bml_matrix_dense_t *bml_zero_matrix_dense_double_complex(
    bml_matrix_dimension_t matrix_dimension,
    bml_distribution_mode_t distrib_mode);

bml_matrix_dense_t *bml_banded_matrix_dense(
    bml_matrix_precision_t matrix_precision,
    int N,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_dense_t *bml_banded_matrix_dense_single_real(
    int N,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_dense_t *bml_banded_matrix_dense_double_real(
    int N,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_dense_t *bml_banded_matrix_dense_single_complex(
    int N,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_dense_t *bml_banded_matrix_dense_double_complex(
    int N,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_dense_t *bml_random_matrix_dense(
    bml_matrix_precision_t matrix_precision,
    int N,
    bml_distribution_mode_t distrib_mode);

bml_matrix_dense_t *bml_random_matrix_dense_single_real(
    int N,
    bml_distribution_mode_t distrib_mode);

bml_matrix_dense_t *bml_random_matrix_dense_double_real(
    int N,
    bml_distribution_mode_t distrib_mode);

bml_matrix_dense_t *bml_random_matrix_dense_single_complex(
    int N,
    bml_distribution_mode_t distrib_mode);

bml_matrix_dense_t *bml_random_matrix_dense_double_complex(
    int N,
    bml_distribution_mode_t distrib_mode);

bml_matrix_dense_t *bml_identity_matrix_dense(
    bml_matrix_precision_t matrix_precision,
    int N,
    bml_distribution_mode_t distrib_mode);

bml_matrix_dense_t *bml_identity_matrix_dense_single_real(
    int N,
    bml_distribution_mode_t distrib_mode);

bml_matrix_dense_t *bml_identity_matrix_dense_double_real(
    int N,
    bml_distribution_mode_t distrib_mode);

bml_matrix_dense_t *bml_identity_matrix_dense_single_complex(
    int N,
    bml_distribution_mode_t distrib_mode);

bml_matrix_dense_t *bml_identity_matrix_dense_double_complex(
    int N,
    bml_distribution_mode_t distrib_mode);

void bml_update_domain_dense(
    bml_matrix_dense_t * A,
    int *localPartMin,
    int *localPartMax,
    int *nnodesInPart);

#endif
