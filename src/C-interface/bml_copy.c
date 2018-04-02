#include "bml_allocate.h"
#include "bml_copy.h"
#include "bml_introspection.h"
#include "bml_parallel.h"
#include "bml_logger.h"
#include "dense/bml_copy_dense.h"
#include "ellpack/bml_copy_ellpack.h"
#include "ellsort/bml_copy_ellsort.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

/** Copy a matrix - result is a new matrix.
 *
 * \ingroup copy_group_C
 *
 * \param A Matrix to copy
 * \return  A Copy of A
 */
bml_matrix_t *
bml_copy_new(
    const bml_matrix_t * A)
{
    bml_matrix_t *B = NULL;

    switch (bml_get_type(A))
    {
        case dense:
            B = bml_copy_dense_new(A);
            break;
        case ellpack:
            B = bml_copy_ellpack_new(A);
            break;
        case ellsort:
            B = bml_copy_ellsort_new(A);
            break;
        default:
            LOG_ERROR("unknown matrix type\n");
            break;
    }
    return B;
}

/** Copy a matrix.
 *
 * \param A Matrix to copy
 * \param B Copy of Matrix A
 */
void
bml_copy(
    const bml_matrix_t * A,
    bml_matrix_t * B)
{
    assert(A != NULL);
    assert(B != NULL);
    LOG_DEBUG("copying matrix\n");
    if (bml_get_type(A) != bml_get_type(B))
    {
        LOG_ERROR("type mismatch\n");
    }
    if (bml_get_N(A) < 0)
    {
        LOG_ERROR("matrix is not initialized\n");
    }
    if (bml_get_N(A) != bml_get_N(B))
    {
        LOG_ERROR("matrix size mismatch\n");
    }
    if (bml_get_M(A) != bml_get_M(B))
    {
        LOG_ERROR("matrix parameter mismatch\n");
    }
    switch (bml_get_type(A))
    {
        case dense:
            bml_copy_dense(A, B);
            break;
        case ellpack:
            bml_copy_ellpack(A, B);
            break;
        case ellsort:
            bml_copy_ellsort(A, B);
            break;
        default:
            LOG_ERROR("unknown matrix type\n");
            break;
    }
}

/** Reorder a matrix in place.
 *
 * \ingroup copy_group_C
 *
 * \param A Matrix to reorder
 * \param perm permutation vector for reordering
 */
void
bml_reorder(
    bml_matrix_t * A,
    int *perm)
{
    switch (bml_get_type(A))
    {
        case dense:
            bml_reorder_dense(A, perm);
            break;
        case ellpack:
            bml_reorder_ellpack(A, perm);
            break;
        case ellsort:
            bml_reorder_ellsort(A, perm);
            break;
        default:
            LOG_ERROR("unknown matrix type\n");
            break;
    }
}

/** Copy a domain.
 *
 * \param A Domain to copy
 * \param B Copy of Domain A
 */
void
bml_copy_domain(
    const bml_domain_t * A,
    bml_domain_t * B)
{
    int nRanks = bml_getNRanks();

    memcpy(B->localRowMin, A->localRowMin, nRanks * sizeof(int));
    memcpy(B->localRowMax, A->localRowMax, nRanks * sizeof(int));
    memcpy(B->localRowExtent, A->localRowExtent, nRanks * sizeof(int));
    memcpy(B->localDispl, A->localDispl, nRanks * sizeof(int));
    memcpy(B->localElements, A->localElements, nRanks * sizeof(int));
}


/** Save current domain for bml matrix.
 *
 * \param A Matrix with domain
 */
void
bml_save_domain(
    bml_matrix_t * A)
{
    switch (bml_get_type(A))
    {
        case dense:
            bml_save_domain_dense(A);
            break;
        case ellpack:
            bml_save_domain_ellpack(A);
            break;
        case ellsort:
            bml_save_domain_ellsort(A);
            break;
        default:
            LOG_ERROR("unknown matrix type\n");
            break;
    }
}

/** Restore to saved domain for bml matrix.
 *
 * \param A Matrix with domain
 */
void
bml_restore_domain(
    bml_matrix_t * A)
{
    switch (bml_get_type(A))
    {
        case dense:
            bml_restore_domain_dense(A);
            break;
        case ellpack:
            bml_restore_domain_ellpack(A);
            break;
        case ellsort:
            bml_restore_domain_ellsort(A);
            break;
        default:
            LOG_ERROR("unknown matrix type\n");
            break;
    }
}
