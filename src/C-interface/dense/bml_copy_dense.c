#include "../bml_allocate.h"
#include "../bml_copy.h"
#include "../bml_logger.h"
#include "../bml_types.h"
#include "bml_allocate_dense.h"
#include "bml_copy_dense.h"
#include "bml_types_dense.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

/** Copy a dense matrix - result in new matrix.
 *
 *  \ingroup copy_group
 *
 *  \param A The matrix to be copied
 *  \return A copy of matrix A.
 */
bml_matrix_dense_t *
bml_copy_dense_new(
    bml_matrix_dense_t * A)
{
    bml_matrix_dense_t *B = NULL;
    assert(A != NULL);
    switch (A->matrix_precision)
    {
        case single_real:
            B = bml_copy_dense_new_single_real(A);
            break;
        case double_real:
            B = bml_copy_dense_new_double_real(A);
            break;
        case single_complex:
            B = bml_copy_dense_new_single_complex(A);
            break;
        case double_complex:
            B = bml_copy_dense_new_double_complex(A);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
    return B;
}

/** Copy a dense matrix.
 *
 *  \ingroup copy_group
 *
 *  \param A The matrix to be copied
 *  \param B Copy of matrix A
 */
void
bml_copy_dense(
    bml_matrix_dense_t * A,
    bml_matrix_dense_t * B)
{
    assert(A != NULL);
    assert(B != NULL);
    switch (A->matrix_precision)
    {
        case single_real:
            bml_copy_dense_single_real(A, B);
            break;
        case double_real:
            bml_copy_dense_double_real(A, B);
            break;
        case single_complex:
            bml_copy_dense_single_complex(A, B);
            break;
        case double_complex:
            bml_copy_dense_double_complex(A, B);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
}

/** Reorder a dense matrix using a permutation vector.
 *
 *  \ingroup copy_group
 *
 *  \param A The matrix to be reordered
 *  \param perm permutaiton matrix for reordering
 */
void
bml_reorder_dense(
    bml_matrix_dense_t * A,
    int *perm)
{
    switch (A->matrix_precision)
    {
        case single_real:
            bml_reorder_dense_single_real(A, perm);
            break;
        case double_real:
            bml_reorder_dense_double_real(A, perm);
            break;
        case single_complex:
            bml_reorder_dense_single_complex(A, perm);
            break;
        case double_complex:
            bml_reorder_dense_double_complex(A, perm);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
}

/** Save the domain for a dense matrix.
 *
 * \ingroup copy_group
 *
 * \param A The matrix with the domain to be saved
 */
void
bml_save_domain_dense(
    bml_matrix_dense_t * A)
{
    bml_copy_domain(A->domain, A->domain2);
}

/** Restore the domain for a dense matrix.
 *
 * \ingroup copy_group
 *
 * \param A The matrix with the domain to be restored
 */
void
bml_restore_domain_dense(
    bml_matrix_dense_t * A)
{
    bml_copy_domain(A->domain2, A->domain);
}
