/** \file */

#ifndef __BML_ALLOCATE_H
#define __BML_ALLOCATE_H

#include "bml_types.h"

#include <stdlib.h>

int bml_allocated(
    bml_matrix_t * A);

void *bml_allocate_memory(
    size_t s);

void *bml_noinit_allocate_memory(
    size_t s);

void bml_free_memory(
    void *ptr);

void bml_free_ptr(
    void **ptr);

void bml_deallocate(
    bml_matrix_t ** A);

void bml_deallocate_domain(
    bml_domain_t * D);

void bml_clear(
    bml_matrix_t * A);

bml_matrix_t *bml_noinit_rectangular_matrix(
    bml_matrix_type_t matrix_type,
    bml_matrix_precision_t matrix_precision,
    bml_matrix_dimension_t matrix_dimension,
    bml_distribution_mode_t distrib_mode);

bml_matrix_t *bml_noinit_matrix(
    bml_matrix_type_t matrix_type,
    bml_matrix_precision_t matrix_precision,
    int N,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_t *bml_zero_matrix(
    bml_matrix_type_t matrix_type,
    bml_matrix_precision_t matrix_precision,
    int N,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_t *bml_random_matrix(
    bml_matrix_type_t matrix_type,
    bml_matrix_precision_t matrix_precision,
    int N,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_t *bml_banded_matrix(
    bml_matrix_type_t matrix_type,
    bml_matrix_precision_t matrix_precision,
    int N,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_t *bml_identity_matrix(
    bml_matrix_type_t matrix_type,
    bml_matrix_precision_t matrix_precision,
    int N,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_domain_t *bml_default_domain(
    int N,
    int M,
    bml_distribution_mode_t distrib_mode);

void bml_update_domain(
    bml_matrix_t * A,
    int *localPartMin,
    int *localPartMax,
    int *nnodesInPart);

#endif
