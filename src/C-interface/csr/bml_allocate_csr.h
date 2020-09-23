#ifndef __BML_ALLOCATE_CSR_H
#define __BML_ALLOCATE_CSR_H

#include "bml_types_csr.h"

/** hash table **/
csr_row_index_hash_t *csr_noinit_table(
    const int alloc_size);

void csr_deallocate_table(
    csr_row_index_hash_t * table);

void csr_reset_table(
    csr_row_index_hash_t * table);

void csr_table_insert(
    csr_row_index_hash_t * table,
    const int key);

void *csr_table_lookup(
    csr_row_index_hash_t * table,
    const int key);

void csr_deallocate_row(
    csr_sparse_row_t * row);

void bml_deallocate_csr(
    bml_matrix_csr_t * A);
/*
void clear_csr_row(
    csr_sparse_row_t * A);
*/

void csr_clear_row_single_real(
    csr_sparse_row_t * A);

void csr_clear_row_double_real(
    csr_sparse_row_t * A);

void csr_clear_row_single_complex(
    csr_sparse_row_t * A);

void csr_clear_row_double_complex(
    csr_sparse_row_t * A);

void bml_clear_csr(
    bml_matrix_csr_t * A);

void bml_clear_csr_single_real(
    bml_matrix_csr_t * A);

void bml_clear_csr_double_real(
    bml_matrix_csr_t * A);

void bml_clear_csr_single_complex(
    bml_matrix_csr_t * A);

void bml_clear_csr_double_complex(
    bml_matrix_csr_t * A);

/*
csr_sparse_row_t *csr_noinit_row(
    const int alloc_size);
*/
csr_sparse_row_t *csr_noinit_row_single_real(
    const int alloc_size);

csr_sparse_row_t *csr_noinit_row_double_real(
    const int alloc_size);

csr_sparse_row_t *csr_noinit_row_single_complex(
    const int alloc_size);

csr_sparse_row_t *csr_noinit_row_double_complex(
    const int alloc_size);

bml_matrix_csr_t *bml_noinit_matrix_csr(
    bml_matrix_precision_t matrix_precision,
    bml_matrix_dimension_t matrix_dimension,
    bml_distribution_mode_t distrib_mode);

bml_matrix_csr_t *bml_noinit_matrix_csr_single_real(
    bml_matrix_dimension_t matrix_dimension,
    bml_distribution_mode_t distrib_mode);

bml_matrix_csr_t *bml_noinit_matrix_csr_double_real(
    bml_matrix_dimension_t matrix_dimension,
    bml_distribution_mode_t distrib_mode);

bml_matrix_csr_t *bml_noinit_matrix_csr_single_complex(
    bml_matrix_dimension_t matrix_dimension,
    bml_distribution_mode_t distrib_mode);

bml_matrix_csr_t *bml_noinit_matrix_csr_double_complex(
    bml_matrix_dimension_t matrix_dimension,
    bml_distribution_mode_t distrib_mode);
/*
csr_sparse_row_t *csr_zero_row (
    const int alloc_size);
*/
csr_sparse_row_t *csr_zero_row_single_real(
    const int alloc_size);

csr_sparse_row_t *csr_zero_row_double_real(
    const int alloc_size);

csr_sparse_row_t *csr_zero_row_single_complex(
    const int alloc_size);

csr_sparse_row_t *csr_zero_row_double_complex(
    const int alloc_size);

bml_matrix_csr_t *bml_zero_matrix_csr(
    bml_matrix_precision_t matrix_precision,
    int N,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_csr_t *bml_zero_matrix_csr_single_real(
    int N,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_csr_t *bml_zero_matrix_csr_double_real(
    int N,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_csr_t *bml_zero_matrix_csr_single_complex(
    int N,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_csr_t *bml_zero_matrix_csr_double_complex(
    int N,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_csr_t *bml_banded_matrix_csr(
    bml_matrix_precision_t matrix_precision,
    int N,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_csr_t *bml_banded_matrix_csr_single_real(
    int N,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_csr_t *bml_banded_matrix_csr_double_real(
    int N,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_csr_t *bml_banded_matrix_csr_single_complex(
    int N,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_csr_t *bml_banded_matrix_csr_double_complex(
    int N,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_csr_t *bml_random_matrix_csr(
    bml_matrix_precision_t matrix_precision,
    int N,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_csr_t *bml_random_matrix_csr_single_real(
    int N,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_csr_t *bml_random_matrix_csr_double_real(
    int N,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_csr_t *bml_random_matrix_csr_single_complex(
    int N,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_csr_t *bml_random_matrix_csr_double_complex(
    int N,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_csr_t *bml_identity_matrix_csr(
    bml_matrix_precision_t matrix_precision,
    int N,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_csr_t *bml_identity_matrix_csr_single_real(
    int N,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_csr_t *bml_identity_matrix_csr_double_real(
    int N,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_csr_t *bml_identity_matrix_csr_single_complex(
    int N,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_csr_t *bml_identity_matrix_csr_double_complex(
    int N,
    int M,
    bml_distribution_mode_t distrib_mode);

void bml_update_domain_csr(
    bml_matrix_csr_t * A,
    int *localPartMin,
    int *localPartMax,
    int *nnodesInPart);

#endif
