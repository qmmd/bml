#include "../../macros.h"
#include "../../typed.h"
#include "../bml_allocate.h"
#include "../bml_introspection.h"
#include "../bml_logger.h"
#include "../bml_types.h"
#include "bml_setters_ellblock.h"
#include "bml_allocate_ellblock.h"
#include "bml_types_ellblock.h"
#include "bml_utilities_ellblock.h"

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <string.h>
#include <assert.h>

/** Set element i,j asuming there's no resetting of any element of A.
 *
 *  \ingroup setters
 *
 *  \param A The matrix which takes row i
 *  \param i The column index
 *  \param j The row index
 *  \param value The element to be added
 *  \WARNING sets an element from scratch
 *  \todo set element new.
 *
 *
 */
void TYPED_FUNC(
    bml_set_element_new_ellblock) (
    bml_matrix_ellblock_t * A,
    int i,
    int j,
    void *_element)
{
    REAL_T *element = _element;
    REAL_T **A_ptr_value = (REAL_T **) A->ptr_value;
    int *A_indexb = A->indexb;
    int *A_nnzb = A->nnzb;
    int *A_bsize = A->bsize;

    //determine block index and index within block
    int ib = 0;
    int jb = 0;
    int ii = i;
    int jj = j;
    while (ii >= A_bsize[ib])
    {
        ii -= A_bsize[ib];
        ib++;
    }
    while (jj >= A_bsize[jb])
    {
        jj -= A_bsize[jb];
        jb++;
    }
    int block_found = 0;
    for (int jp = 0; jp < A_nnzb[ib]; jp++)
    {
        int ind = ROWMAJOR(ib, jp, A->NB, A->MB);
        if (A_indexb[ind] == jb)
        {
            REAL_T *A_value = A_ptr_value[ind];
            assert(A_value != NULL);
            A_value[ROWMAJOR(ii, jj, A_bsize[ib], A_bsize[jb])]
                = *((REAL_T *) element);
            block_found = 1;
        }
    }
    if (block_found == 0)
    {
        if (A_nnzb[ib] == A->MB)
        {
            LOG_ERROR("Number of non-zeroes per row > MB, Increase MB\n");
        }
        int ind = ROWMAJOR(ib, A_nnzb[ib], A->NB, A->MB);
        assert(A_ptr_value[ind] == NULL);

        // printf("allocate new block ib=%d, ind=%d\n", ib, ind);
        int nelements = A_bsize[ib] * A_bsize[jb];
        A_ptr_value[ind] =
            TYPED_FUNC(bml_allocate_block_ellblock) (A, ib, nelements);
        A_indexb[ind] = jb;

        REAL_T *A_value = A_ptr_value[ind];
        assert(A_value != NULL);
        memset(A_value, 0, nelements * sizeof(REAL_T));
        A_value[ROWMAJOR(ii, jj, A_bsize[ib], A_bsize[jb])]
            = *((REAL_T *) element);
        A_nnzb[ib]++;
    }
}


/** Set element i,j of matrix A.
 *
 *  \ingroup setters
 *
 *  \param A The matrix which takes row i
 *  \param i The column index
 *  \param j The row index
 *  \param value The element to be set
 *  \WARNING sets an element from scratch
 *  \todo set element new.
 *
 *
 */
void TYPED_FUNC(
    bml_set_element_ellblock) (
    bml_matrix_ellblock_t * A,
    int i,
    int j,
    void *element)
{
    TYPED_FUNC(bml_set_element_new_ellblock) (A, i, j, element);
}

/** Set row i of matrix A.
 *
 *  \ingroup setters
 *
 *  \param A The matrix which takes row i
 *  \param i The index of the row to be set
 *  \param row The row to be set
 *  \param threshold The threshold value to be set
 *
 */
void TYPED_FUNC(
    bml_set_row_ellblock) (
    bml_matrix_ellblock_t * A,
    int i,
    void *_row,
    double threshold)
{
    REAL_T *row = _row;
    int A_N = A->N;

    REAL_T **A_ptr_value = (REAL_T **) A->ptr_value;
    int *A_indexb = A->indexb;
    int *A_nnzb = A->nnzb;
    int *A_bsize = A->bsize;

    //determine block index and index within block
    int ib = 0;
    int ii = i;
    while (ii >= A_bsize[ib])
    {
        ii -= A_bsize[ib];
        ib++;
    }

    //loop over elements to insert
    for (int k = 0; k < A_N; k++)
    {
        // printf("k=%d\n", k);
        if (ABS(row[k]) > threshold)
        {
            //determine column block index and index within block
            int jb = 0;
            int jj = k;
            while (jj >= A_bsize[jb])
            {
                jj -= A_bsize[jb];
                jb++;
            }
            int block_found = 0;
            for (int jp = 0; jp < A_nnzb[ib]; jp++)
            {
                int ind = ROWMAJOR(ib, jp, A->NB, A->MB);
                if (A_indexb[ind] == jb)
                {
                    REAL_T *A_value = A_ptr_value[ind];
                    A_value[ROWMAJOR(ii, jj, A_bsize[ib], A_bsize[jb])]
                        = row[k];
                    block_found = 1;
                }
            }
            if (block_found == 0)
            {
                int ind = ROWMAJOR(ib, A_nnzb[ib], A->NB, A->MB);
                int nelements = A_bsize[ib] * A_bsize[jb];
                A_ptr_value[ind]
                    =
                    TYPED_FUNC(bml_allocate_block_ellblock) (A, ib,
                                                             nelements);
                REAL_T *A_value = A_ptr_value[ind];
                memset(A_value, 0, nelements * sizeof(REAL_T));
                assert(A_value != NULL);
                A_value[ROWMAJOR(ii, jj, A_bsize[ib], A_bsize[jb])] = row[k];
                A_nnzb[ib]++;
                A_indexb[ind] = jb;
            }
        }
    }

}

/** Set diagonal of matrix A.
 *
 *  \ingroup setters
 *
 *  \param A The matrix which takes diag
 *  \param diag The diagonal to be set
 *  \param threshold The threshold value to be used
 */
void TYPED_FUNC(
    bml_set_diagonal_ellblock) (
    bml_matrix_ellblock_t * A,
    void *_diagonal,
    double threshold)
{
    REAL_T *diagonal = _diagonal;
    int A_NB = A->NB;
    int A_MB = A->MB;

    REAL_T **A_ptr_value = (REAL_T **) A->ptr_value;
    int *A_indexb = A->indexb;
    int *A_nnzb = A->nnzb;
    int *A_bsize = A->bsize;

    int offset = 0;
    for (int ib = 0; ib < A_NB; ib++)
    {
        int block_found = 0;
        for (int jp = 0; jp < A_nnzb[ib]; jp++)
        {
            int ind = ROWMAJOR(ib, jp, A_NB, A_MB);
            int jb = A_indexb[ind];
            if (jb == ib)
            {
                REAL_T *A_value = A_ptr_value[ind];
                double normdiag = TYPED_FUNC(bml_norm_inf)
                    (&diagonal[offset], A_bsize[ib], 1, 1);

                if (normdiag > threshold)
                {
                    for (int ii = 0; ii < A_bsize[ib]; ii++)
                    {
                        A_value[ROWMAJOR(ii, ii, A_bsize[ib], A_bsize[ib])]
                            = diagonal[offset + ii];
                    }
                }
                else
                {
                    for (int ii = 0; ii < A_bsize[ib]; ii++)
                    {
                        A_value[ROWMAJOR(ii, ii, A_bsize[ib], A_bsize[ib])] =
                            0.0;
                    }
                }
                block_found = 1;
            }
        }

        /* If there is no diagonal block then
         */
        if (block_found == 0)
        {
            double normdiag = TYPED_FUNC(bml_norm_inf)
                (&diagonal[offset], A_bsize[ib], 1, 1);
            if (normdiag > threshold)
            {
                int ind = ROWMAJOR(ib, A_nnzb[ib], A->NB, A->MB);
                int nelemnts = A_bsize[ib] * A_bsize[ib];
                A_ptr_value[ind] =
                    TYPED_FUNC(bml_allocate_block_ellblock) (A, ib, nelemnts);
                A_indexb[ind] = ib;

                REAL_T *A_value = A_ptr_value[ind];
                memset(A_value, 0, nelemnts * sizeof(REAL_T));
                for (int ii = 0; ii < A_bsize[ib]; ii++)
                {
                    A_value[ROWMAJOR(ii, ii, A_bsize[ib], A_bsize[ib])]
                        = diagonal[offset + ii];
                }
                A_nnzb[ib]++;
            }
        }
        offset += A_bsize[ib];
    }
}

/*
 * This function assumes the block structure of the matrix
 * has already been set and the block ib, jb already exists
 */
void TYPED_FUNC(
    bml_set_block_ellblock) (
    bml_matrix_ellblock_t * A,
    int ib,
    int jb,
    void *_elements)
{
    assert(ib < A->NB);
    assert(jb < A->NB);

    REAL_T *elements = _elements;

    int data_copied = 0;
    for (int jp = 0; jp < A->nnzb[ib]; jp++)
    {
        int ind = ROWMAJOR(ib, jp, A->NB, A->MB);
        if (A->indexb[ind] == jb)
        {
            int n2 = A->bsize[ib] * A->bsize[jb];
            REAL_T *A_value = A->ptr_value[ind];
            assert(A_value != NULL);
            memcpy(A_value, elements, n2 * sizeof(REAL_T));
            data_copied = 1;
            break;
        }
    }

    //make sure block was found and data was copied
    assert(data_copied == 1);
}
