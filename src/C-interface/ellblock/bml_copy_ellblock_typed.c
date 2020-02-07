#include "../../macros.h"
#include "../../typed.h"
#include "../bml_allocate.h"
#include "../bml_copy.h"
#include "../bml_types.h"
#include "bml_allocate_ellblock.h"
#include "bml_copy_ellblock.h"
#include "bml_types_ellblock.h"

#include <assert.h>
#include <complex.h>
#include <stdlib.h>
#include <string.h>

/** Copy an ellblock matrix - result is a new matrix.
 *
 *  \ingroup copy_group
 *
 *  \param A The matrix to be copied
 *  \return A copy of matrix A.
 */
bml_matrix_ellblock_t
    * TYPED_FUNC(bml_copy_ellblock_new) (bml_matrix_ellblock_t * A)
{
    bml_matrix_ellblock_t *B =
        TYPED_FUNC(bml_block_matrix_ellblock) (A->NB, A->MB, A->bsize,
                                               A->distribution_mode);

    int NB = A->NB;
    int MB = A->MB;

    int *A_indexb = A->indexb;
    REAL_T **A_ptr_value = (REAL_T **) A->ptr_value;

    int *B_indexb = B->indexb;
    REAL_T **B_ptr_value = (REAL_T **) B->ptr_value;

    //    memcpy(B->index, A->index, sizeof(int) * A->N * A->M);
    memcpy(B->nnzb, A->nnzb, sizeof(int) * A->NB);
    //    memcpy(B->value, A->value, sizeof(REAL_T) * A->N * A->M);
    memcpy(B->bsize, A->bsize, sizeof(int) * A->NB);

    memcpy(B_indexb, A_indexb, NB * MB * sizeof(int));

    for (int ib = 0; ib < NB; ib++)
    {
        for (int jp = 0; jp < A->nnzb[ib]; jp++)
        {
            int ind = ROWMAJOR(ib, jp, NB, MB);
            int jb = B_indexb[ind];
            B_ptr_value[ind]
                = malloc(A->bsize[ib] * A->bsize[jb] * sizeof(REAL_T));

            memcpy(B_ptr_value[ind], A_ptr_value[ind],
                   A->bsize[ib] * A->bsize[jb] * sizeof(REAL_T));
        }
    }

    return B;
}

/** Copy an ellblock matrix.
 *
 *  \ingroup copy_group
 *
 *  \param A The matrix to be copied
 *  \param B Copy of matrix A
 */
void TYPED_FUNC(
    bml_copy_ellblock) (
    bml_matrix_ellblock_t * A,
    bml_matrix_ellblock_t * B)
{
    assert(A->NB == B->NB);

    int NB = A->NB;
    int MB = A->MB;

    REAL_T **A_ptr_value = (REAL_T **) A->ptr_value;
    REAL_T **B_ptr_value = (REAL_T **) B->ptr_value;

    memcpy(B->nnzb, A->nnzb, sizeof(int) * A->NB);

    int *A_indexb = A->indexb;
    int *B_indexb = B->indexb;
    memcpy(B_indexb, A_indexb, NB * MB * sizeof(int));

    for (int ib = 0; ib < NB; ib++)
    {
        for (int jp = 0; jp < A->nnzb[ib]; jp++)
        {
            int ind = ROWMAJOR(ib, jp, NB, MB);
            assert(A_ptr_value[ind] != NULL);
            int jb = B_indexb[ind];
            int nelements = A->bsize[ib] * A->bsize[jb];
            if (B_ptr_value[ind] == NULL)
                B_ptr_value[ind]
                    = bml_noinit_allocate_memory(nelements * sizeof(REAL_T));
            memcpy(B_ptr_value[ind], A_ptr_value[ind],
                   nelements * sizeof(REAL_T));
        }
    }
}

/** Reorder an ellblock matrix.
 *
 *  \ingroup copy_group
 *
 *  \param A The matrix to be reordered
 *  \param B The permutation vector
 */
void TYPED_FUNC(
    bml_reorder_ellblock) (
    bml_matrix_ellblock_t * A,
    int *perm)
{
    int NB = A->NB;
    int MB = A->MB;

    int *A_indexb = A->indexb;
    int *A_nnzb = A->nnzb;
    REAL_T **A_ptr_value = (REAL_T **) A->ptr_value;
    int *A_bsize = A->bsize;

    bml_matrix_ellblock_t *B = bml_copy_new(A);
    int *B_indexb = B->indexb;
    int *B_nnzb = B->nnzb;
    REAL_T **B_ptr_value = (REAL_T **) B->ptr_value;
    int *B_bsize = B->bsize;

    for (int i = 0; i < NB; i++)
    {
        A_bsize[i] = B_bsize[perm[i]];
    }

    // Reorder rows - need to copy
    for (int ib = 0; ib < NB; ib++)
    {
        memcpy(&A_indexb[ROWMAJOR(perm[ib], 0, NB, MB)],
               &B_indexb[ROWMAJOR(ib, 0, NB, MB)], MB * sizeof(int));
        int count = 0;
        for (int jp = 0; jp < MB; jp++)
        {
            int ind = ROWMAJOR(ib, jp, NB, MB);
            int jb = B_indexb[ind];
            count += B_bsize[jb];
        }
        memcpy(A_ptr_value[ROWMAJOR(perm[ib], 0, NB, MB)],
               B_ptr_value[ROWMAJOR(ib, 0, NB, MB)],
               B_bsize[ib] * count * sizeof(REAL_T));
        A_nnzb[perm[ib]] = B_nnzb[ib];
    }

    bml_deallocate_ellblock(B);

    // Reorder elements in each row - just change index
    for (int ib = 0; ib < NB; ib++)
    {
        for (int jp = 0; jp < A_nnzb[ib]; jp++)
        {
            A_indexb[ROWMAJOR(ib, jp, NB, MB)] =
                perm[A_indexb[ROWMAJOR(ib, jp, NB, MB)]];
        }
    }
}
