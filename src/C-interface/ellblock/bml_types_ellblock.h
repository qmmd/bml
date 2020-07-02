#ifndef __BML_TYPES_ELLBLOCK_H
#define __BML_TYPES_ELLBLOCK_H

#include "../bml_types.h"

/** BLOCK ELLPACK matrix type. */
struct bml_matrix_ellblock_t
{
    /** The matrix type identifier. */
    bml_matrix_type_t matrix_type;
    /** The real precision. */
    bml_matrix_precision_t matrix_precision;
    /** The distribution mode. **/
    bml_distribution_mode_t distribution_mode;
    /** The number of rows. */
    int N;
    /** The number of columns per row. */
    int M;
    /** The number of block rows. */
    int NB;
    /** The max. number of blocks per row. */
    int MB;
#ifdef BML_ELLBLOCK_USE_MEMPOOL
    /** Storage for matrix elements/blocks. */
    void *memory_pool;
    /** offsets for storage of row blocks */
    int *memory_pool_offsets;
    /** Ptr to last stored block in row block */
    void **memory_pool_ptr;
#endif
    /** Pointers to blocks of values */
    void **ptr_value;
    /** The block index matrix. */
    int *indexb;
    /** The vector of non-zeros blocks per row */
    int *nnzb;
    /** The dimensions of the blocks */
    int *bsize;
    /** The domain decomposition when running in parallel. */
    bml_domain_t *domain;
    /** A copy of the domain decomposition. */
    bml_domain_t *domain2;
};
typedef struct bml_matrix_ellblock_t bml_matrix_ellblock_t;

#define BMAXSIZE 25

#endif
