/** \file */

#ifndef __BML_ALLOCATE_H
#define __BML_ALLOCATE_H

#include "bml_types.h"

#include <stdlib.h>

void *bml_allocate_memory(
    const size_t s);

void *bml_noinit_allocate_memory(
    const size_t s);

void bml_free_memory(
    void *ptr);

void bml_free_ptr(
    double **ptr);

void bml_deallocate(
    bml_matrix_t ** A);

void bml_clear(
    bml_matrix_t * A);

int bml_allocated(
    const bml_matrix_t * A);

bml_matrix_t *bml_noinit_matrix(
    const bml_matrix_type_t matrix_type,
    const bml_matrix_precision_t matrix_precision,
    const int N,
    const int M,
    const bml_distribution_mode_t distrib_mode);

bml_matrix_t *bml_zero_matrix(
    const bml_matrix_type_t matrix_type,
    const bml_matrix_precision_t matrix_precision,
    const int N,
    const int M,
    const bml_distribution_mode_t distrib_mode);

bml_matrix_t *bml_banded_matrix(
    const bml_matrix_type_t matrix_type,
    const bml_matrix_precision_t matrix_precision,
    const int N,
    const int M,
    const bml_distribution_mode_t distrib_mode);

bml_matrix_t *bml_random_matrix(
    const bml_matrix_type_t matrix_type,
    const bml_matrix_precision_t matrix_precision,
    const int N,
    const int M,
    const bml_distribution_mode_t distrib_mode);

bml_matrix_t *bml_identity_matrix(
    const bml_matrix_type_t matrix_type,
    const bml_matrix_precision_t matrix_precision,
    const int N,
    const int M,
    const bml_distribution_mode_t distrib_mode);

void bml_deallocate_domain(
    bml_domain_t * D);

bml_domain_t *bml_default_domain(
    const int N,
    const int M,
    const bml_distribution_mode_t distrib_mode);

void bml_update_domain(
    bml_matrix_t * A,
    int *localPartMin,
    int *localPartMax,
    int *nnodesInPart);

#endif
