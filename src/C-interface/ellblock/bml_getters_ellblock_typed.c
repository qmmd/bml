#include "bml_getters_ellblock.h"
#include "../bml_introspection.h"
#include "../bml_logger.h"
#include "bml_types_ellblock.h"
#include "../../macros.h"
#include "../../typed.h"

#include <complex.h>
#include <stdlib.h>
#include <assert.h>

/** Return a single matrix element.
 *
 * \param A The bml matrix
 * \param i The row index
 * \param j The column index
 * \return The matrix element
 */
void *TYPED_FUNC(
    bml_get_ellblock) (
    bml_matrix_ellblock_t * A,
    int i,
    int j)
{
    static REAL_T MINUS_ONE = -1;
    static REAL_T ZERO = 0;
    REAL_T **A_ptr_value = (REAL_T **) A->ptr_value;
    int *bsize = A->bsize;

    if (i < 0 || i >= A->N)
    {
        LOG_ERROR("row index out of bounds\n");
        return &MINUS_ONE;
    }
    if (j < 0 || j >= A->N)
    {
        LOG_ERROR("column index out of bounds\n");
        return &MINUS_ONE;
    }
    //determine block index and index within block
    int ib = 0;
    int jb = 0;
    int ii = i;
    int jj = j;
    while (ii >= bsize[ib])
    {
        ii -= bsize[ib];
        ib++;
    }
    while (jj >= bsize[jb])
    {
        jj -= bsize[jb];
        jb++;
    }
    for (int jp = 0; jp < A->nnzb[ib]; jp++)
    {
        int ind = ROWMAJOR(ib, jp, A->NB, A->MB);
        if (A->indexb[ind] == jb)
        {
            REAL_T *A_value = A_ptr_value[ind];
            return &A_value[ROWMAJOR(ii, jj, bsize[ib], bsize[jb])];
        }
    }
    return &ZERO;
}

/** Get row i of matrix A.
 *
 *  \ingroup getters
 *
 *  \param A The matrix which takes row i
 *  \param i The index of the row to get
 *  \param row Array to copy the row
 *
 */
void *TYPED_FUNC(
    bml_get_row_ellblock) (
    bml_matrix_ellblock_t * A,
    int i)
{
    int A_NB = A->NB;
    int A_MB = A->MB;
    REAL_T **A_ptr_value = (REAL_T **) A->ptr_value;
    int *A_indexb = A->indexb;
    int *A_nnzb = A->nnzb;
    int *bsize = A->bsize;

    REAL_T *row = calloc(A->N, sizeof(REAL_T));
    for (int ii = 0; ii < A->N; ii++)
    {
        row[ii] = 0.0;
    }

    //determine row block index and row index within block
    int ib = 0;
    int ii = i;
    while (ii >= bsize[ib])
    {
        ii -= bsize[ib];
        ib++;
    }

    int *offset = malloc(A_NB * sizeof(int));
    offset[0] = 0;
    for (int jb = 1; jb < A_NB; jb++)
        offset[jb] = offset[jb - 1] + bsize[jb - 1];

    for (int jp = 0; jp < A_nnzb[ib]; jp++)
    {
        int ind = ROWMAJOR(ib, jp, A_NB, A_MB);
        int jb = A_indexb[ind];
        if (jb >= 0)
        {
            REAL_T *A_value = A_ptr_value[ind];
            for (int jj = 0; jj < bsize[jb]; jj++)
            {
                int ll = offset[jb] + jj;
                row[ll] = A_value[ROWMAJOR(ii, jj, bsize[ib], bsize[jb])];
            }
        }
    }
    free(offset);

    return row;
}

/** Get the diagonal of matrix A.
 *
 *  \ingroup getters
 *
 *  \param A The matrix which takes row i
 *  \param Diagonal Array to copy the diagonal
 *
 */
void *TYPED_FUNC(
    bml_get_diagonal_ellblock) (
    bml_matrix_ellblock_t * A)
{
    int NB = A->NB;
    int MB = A->MB;
    REAL_T **A_ptr_value = (REAL_T **) A->ptr_value;
    int *indexb = A->indexb;
    int *nnzb = A->nnzb;
    int *bsize = A->bsize;
    REAL_T *diagonal = calloc(A->N, sizeof(REAL_T));

    int offset = 0;
    for (int ib = 0; ib < NB; ib++)
    {
        for (int jp = 0; jp < nnzb[ib]; jp++)
        {
            int ind = ROWMAJOR(ib, jp, NB, MB);
            if (indexb[ind] == ib)
            {
                REAL_T *A_value = A_ptr_value[ind];
                assert(A_value != NULL);
                for (int ii = 0; ii < bsize[ib]; ii++)
                    diagonal[offset + ii]
                        = A_value[ROWMAJOR(ii, ii, bsize[ib], bsize[ib])];
            }
        }
        offset += bsize[ib];
    }

    return diagonal;
}
