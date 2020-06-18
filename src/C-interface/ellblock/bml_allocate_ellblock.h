#ifndef __BML_ALLOCATE_ELLBLOCK_H
#define __BML_ALLOCATE_ELLBLOCK_H

#include "bml_types_ellblock.h"

void bml_deallocate_ellblock(
    bml_matrix_ellblock_t * A);

void bml_clear_ellblock(
    bml_matrix_ellblock_t * A);

void bml_clear_ellblock_single_real(
    bml_matrix_ellblock_t * A);

void bml_clear_ellblock_double_real(
    bml_matrix_ellblock_t * A);

void bml_clear_ellblock_single_complex(
    bml_matrix_ellblock_t * A);

void bml_clear_ellblock_double_complex(
    bml_matrix_ellblock_t * A);

bml_matrix_ellblock_t *bml_noinit_matrix_ellblock(
    bml_matrix_precision_t matrix_precision,
    bml_matrix_dimension_t matrix_dimension,
    bml_distribution_mode_t distrib_mode);

bml_matrix_ellblock_t
    * bml_noinit_matrix_ellblock_single_real(bml_matrix_dimension_t
                                             matrix_dimension,
                                             bml_distribution_mode_t
                                             distrib_mode);

bml_matrix_ellblock_t
    * bml_noinit_matrix_ellblock_double_real(bml_matrix_dimension_t
                                             matrix_dimension,
                                             bml_distribution_mode_t
                                             distrib_mode);

bml_matrix_ellblock_t
    * bml_noinit_matrix_ellblock_single_complex(bml_matrix_dimension_t
                                                matrix_dimension,
                                                bml_distribution_mode_t
                                                distrib_mode);

bml_matrix_ellblock_t
    * bml_noinit_matrix_ellblock_double_complex(bml_matrix_dimension_t
                                                matrix_dimension,
                                                bml_distribution_mode_t
                                                distrib_mode);

bml_matrix_ellblock_t *bml_zero_matrix_ellblock(
    bml_matrix_precision_t matrix_precision,
    int N,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_ellblock_t *bml_zero_matrix_ellblock_single_real(
    int N,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_ellblock_t *bml_zero_matrix_ellblock_double_real(
    int N,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_ellblock_t *bml_zero_matrix_ellblock_single_complex(
    int N,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_ellblock_t *bml_zero_matrix_ellblock_double_complex(
    int N,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_ellblock_t *bml_banded_matrix_ellblock(
    bml_matrix_precision_t matrix_precision,
    int N,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_ellblock_t *bml_banded_matrix_ellblock_single_real(
    int N,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_ellblock_t *bml_banded_matrix_ellblock_double_real(
    int N,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_ellblock_t *bml_banded_matrix_ellblock_single_complex(
    int N,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_ellblock_t *bml_banded_matrix_ellblock_double_complex(
    int N,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_ellblock_t *bml_random_matrix_ellblock(
    bml_matrix_precision_t matrix_precision,
    int N,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_ellblock_t *bml_random_matrix_ellblock_single_real(
    int N,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_ellblock_t *bml_random_matrix_ellblock_double_real(
    int N,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_ellblock_t *bml_random_matrix_ellblock_single_complex(
    int N,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_ellblock_t *bml_random_matrix_ellblock_double_complex(
    int N,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_ellblock_t *bml_identity_matrix_ellblock(
    bml_matrix_precision_t matrix_precision,
    int N,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_ellblock_t *bml_identity_matrix_ellblock_single_real(
    int N,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_ellblock_t *bml_identity_matrix_ellblock_double_real(
    int N,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_ellblock_t *bml_identity_matrix_ellblock_single_complex(
    int N,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_ellblock_t *bml_identity_matrix_ellblock_double_complex(
    int N,
    int M,
    bml_distribution_mode_t distrib_mode);

bml_matrix_ellblock_t *bml_block_matrix_ellblock(
    bml_matrix_precision_t matrix_precision,
    int NB,
    int MB,
    int M,
    int *bsizes,
    bml_distribution_mode_t distrib_mode);

bml_matrix_ellblock_t *bml_block_matrix_ellblock_single_real(
    int NB,
    int MB,
    int M,
    int *bsizes,
    bml_distribution_mode_t distrib_mode);

bml_matrix_ellblock_t *bml_block_matrix_ellblock_double_real(
    int NB,
    int MB,
    int M,
    int *bsizes,
    bml_distribution_mode_t distrib_mode);

bml_matrix_ellblock_t *bml_block_matrix_ellblock_single_complex(
    int NB,
    int MB,
    int M,
    int *bsizes,
    bml_distribution_mode_t distrib_mode);

bml_matrix_ellblock_t *bml_block_matrix_ellblock_double_complex(
    int NB,
    int MB,
    int M,
    int *bsizes,
    bml_distribution_mode_t distrib_mode);

void bml_update_domain_ellblock(
    bml_matrix_ellblock_t * A,
    int *localPartMin,
    int *localPartMax,
    int *nnodesInPart);

void bml_set_block_sizes(
    int *bsize,
    int NB,
    int MB);
int *bml_get_block_sizes(
    int N,
    int M);
int bml_get_mb(
    );
int bml_get_nb(
    );
int count_nelements(
    bml_matrix_ellblock_t * A);

void *bml_allocate_block_ellblock_single_real(
    bml_matrix_ellblock_t * A,
    const int ib,
    const int nelements);
void *bml_allocate_block_ellblock_double_real(
    bml_matrix_ellblock_t * A,
    const int ib,
    const int nelements);
void *bml_allocate_block_ellblock_single_complex(
    bml_matrix_ellblock_t * A,
    const int ib,
    const int nelements);
void *bml_allocate_block_ellblock_double_complex(
    bml_matrix_ellblock_t * A,
    const int ib,
    const int nelements);

void bml_free_block_ellblock_single_real(
    bml_matrix_ellblock_t * A,
    const int ib,
    const int jb);
void bml_free_block_ellblock_double_real(
    bml_matrix_ellblock_t * A,
    const int ib,
    const int jb);
void bml_free_block_ellblock_single_complex(
    bml_matrix_ellblock_t * A,
    const int ib,
    const int jb);
void bml_free_block_ellblock_double_complex(
    bml_matrix_ellblock_t * A,
    const int ib,
    const int jb);
#endif
