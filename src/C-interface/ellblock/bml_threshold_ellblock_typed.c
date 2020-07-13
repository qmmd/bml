#include "../../macros.h"
#include "../../typed.h"
#include "../bml_allocate.h"
#include "../bml_parallel.h"
#include "../bml_threshold.h"
#include "../bml_types.h"
#include "bml_allocate_ellblock.h"
#include "bml_threshold_ellblock.h"
#include "bml_types_ellblock.h"
#include "bml_utilities_ellblock.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/** Threshold a matrix.
 *
 *  \ingroup threshold_group
 *
 *  \param A The matrix to be thresholded
 *  \param threshold Threshold value
 *  \return the thresholded A
 */
bml_matrix_ellblock_t
    * TYPED_FUNC(bml_threshold_new_ellblock) (bml_matrix_ellblock_t * A,
                                              double threshold)
{
    int NB = A->NB;
    int MB = A->MB;
    int *bsize = A->bsize;

    bml_matrix_ellblock_t *B =
        TYPED_FUNC(bml_block_matrix_ellblock) (NB, MB, A->M, bsize,
                                               A->distribution_mode);

    REAL_T **A_ptr_value = (REAL_T **) A->ptr_value;
    int *A_indexb = A->indexb;
    int *A_nnzb = A->nnzb;

    REAL_T **B_ptr_value = (REAL_T **) B->ptr_value;
    int *B_indexb = B->indexb;
    int *B_nnzb = B->nnzb;

    for (int ib = 0; ib < NB; ib++)
    {
        for (int jp = 0; jp < A_nnzb[ib]; jp++)
        {
            int ind = ROWMAJOR(ib, jp, NB, MB);
            int jb = A_indexb[ind];
            REAL_T *A_value = A_ptr_value[ind];
            REAL_T normA = TYPED_FUNC(bml_norm_inf)
                (A_value, bsize[ib], bsize[jb], bsize[jb]);

            if (is_above_threshold(normA, threshold))
            {
                int nelements = bsize[ib] * bsize[jb];
                int indB = ROWMAJOR(ib, B_nnzb[ib], NB, MB);
                B_ptr_value[indB]
                    =
                    TYPED_FUNC(bml_allocate_block_ellblock) (B, ib,
                                                             nelements);
                REAL_T *B_value = B_ptr_value[indB];
                memcpy(B_value, A_value, nelements * sizeof(REAL_T));
                for (int ii = 0; ii < bsize[ib]; ii++)
                    for (int jj = 0; jj < bsize[jb]; jj++)
                    {
                        int index = ROWMAJOR(ii, jj, bsize[ib], bsize[jb]);
                        if (!is_above_threshold(B_value[index], threshold))
                        {
                            B_value[index] = 0.;
                        }
                    }
                B_indexb[indB] = jb;
                B_nnzb[ib]++;
            }
        }
    }

    return B;
}

/** Threshold a matrix in place.
 *
 *  \ingroup threshold_group
 *
 *  \param A The matrix to be thresholded
 *  \param threshold Threshold value
 *  \return the thresholded A
 */
void TYPED_FUNC(
    bml_threshold_ellblock) (
    bml_matrix_ellblock_t * A,
    double threshold)
{
    int NB = A->NB;
    int MB = A->MB;

    REAL_T **A_ptr_value = (REAL_T **) A->ptr_value;
    int *A_indexb = A->indexb;
    int *A_nnzb = A->nnzb;
    int *bsize = A->bsize;

    int rlen;

    for (int ib = 0; ib < NB; ib++)
    {
        for (int jp = 0; jp < A_nnzb[ib]; jp++)
        {
            int ind = ROWMAJOR(ib, jp, NB, MB);
            int jb = A_indexb[ind];
            REAL_T *A_value = A_ptr_value[ind];
            REAL_T normA = TYPED_FUNC(bml_norm_inf)
                (A_value, bsize[ib], bsize[jb], bsize[jb]);

            if (is_above_threshold(normA, threshold))
            {
                //apply thresholding within block
                REAL_T *B_value = A_ptr_value[ind];
                for (int ii = 0; ii < bsize[ib]; ii++)
                    for (int jj = 0; jj < bsize[jb]; jj++)
                    {
                        int index = ROWMAJOR(ii, jj, bsize[ib], bsize[jb]);
                        if (!is_above_threshold(B_value[index], threshold))
                        {
                            B_value[index] = 0.;
                        }
                    }
            }
            else
            {
                //printf("remove block %d %d\n",ib,jb);
                TYPED_FUNC(bml_free_block_ellblock) (A, ib, jb);
                jp--;
            }
        }
    }
}
