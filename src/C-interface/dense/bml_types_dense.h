#ifndef __BML_TYPES_DENSE_H
#define __BML_TYPES_DENSE_H

#include "../bml_types.h"

#ifdef BML_USE_MAGMA
#include "magma_auxiliary.h"
#endif

/** Dense matrix type. */
struct bml_matrix_dense_t
{
    /** The matrix type identifier. */
    bml_matrix_type_t matrix_type;
    /** The real precision. */
    bml_matrix_precision_t matrix_precision;
    /** The distribution mode. **/
    bml_distribution_mode_t distribution_mode;
    /** The number of rows/columns. */
    int N;
    /** The dense matrix. */
    void *matrix;
    /** The leading dimension of the array matrix. */
    int ld;
    /** The domain decomposition when running in parallel. */
    bml_domain_t *domain;
    /** A copy of the domain decomposition. */
    bml_domain_t *domain2;
#ifdef BML_USE_MAGMA
    magma_queue_t queue;
#endif
};
typedef struct bml_matrix_dense_t bml_matrix_dense_t;

#endif
