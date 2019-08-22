#include "../../macros.h"
#include "../../typed.h"
#include "../bml_allocate.h"
#include "../bml_parallel.h"
#include "../bml_transpose.h"
#include "../bml_types.h"
#include "bml_allocate_ellblock.h"
#include "bml_transpose_ellblock.h"
#include "bml_types_ellblock.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/** Transpose a matrix.
 *
 *  \ingroup transpose_group
 *
 *  \param A The matrix to be transposed
 *  \return the transposed A
 */
bml_matrix_ellblock_t *TYPED_FUNC(
    bml_transpose_new_ellblock) (
    bml_matrix_ellblock_t * A)
{
    int NB = A->NB;
    int MB = A->MB;
    int *bsize = A->bsize;
    bml_distribution_mode_t distmode = A->distribution_mode;

    bml_matrix_ellblock_t *B = TYPED_FUNC(bml_block_matrix_ellblock) (NB, MB,
                                                                      bsize,
                                                                      distmode);

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
            int indA = ROWMAJOR(ib, jp, NB, MB);
            int jb = A_indexb[indA];
            int indB = ROWMAJOR(jb, B_nnzb[jb], NB, MB);
            int nelements = bsize[ib] * bsize[jb];
            //add a new block in B
            B_indexb[indB] = ib;
            B_ptr_value[indB] =
                bml_noinit_allocate_memory(nelements * sizeof(REAL_T));
            B_nnzb[jb]++;
            //transpose block
            REAL_T *B_value = B_ptr_value[indB];
            REAL_T *A_value = A_ptr_value[indA];
            for (int ii = 0; ii < bsize[ib]; ii++)
                for (int jj = 0; jj < bsize[jb]; jj++)
                    B_value[ROWMAJOR(jj, ii, bsize[jb], bsize[ib])] =
                        A_value[ROWMAJOR(ii, jj, bsize[ib], bsize[jb])];
        }
    }

    return B;
}


/** Transpose a matrix in place.
 *
 *  \ingroup transpose_group
 *
 *  \param A The matrix to be transposeed
 *  \return the transposed A
 */
void TYPED_FUNC(
    bml_transpose_ellblock) (
    bml_matrix_ellblock_t * A)
{
    int NB = A->NB;
    int MB = A->MB;

    REAL_T **A_ptr_value = (REAL_T **) A->ptr_value;
    int *A_indexb = A->indexb;
    int *A_nnzb = A->nnzb;
    int *bsize = A->bsize;

    for (int ib = 0; ib < NB; ib++)
    {
        for (int jp = A_nnzb[ib] - 1; jp >= 0; jp--)
        {
            int indl = ROWMAJOR(ib, jp, NB, MB);
            int jb = A_indexb[indl];
            if (jb >= ib)
            {
                int exchangeDone = 0;
                for (int kp = 0; kp < A_nnzb[jb]; kp++)
                {
                    int indr = ROWMAJOR(jb, kp, NB, MB);
                    if (A_indexb[indr] == ib)
                    {
                        REAL_T *A_value_l = A_ptr_value[indl];
                        REAL_T *A_value_r = A_ptr_value[indr];
                        if (ib == jb)
                        {
                            for (int ii = 0; ii < bsize[ib]; ii++)
                                for (int jj = 0; jj < ii; jj++)
                                {
                                    int il = ROWMAJOR(ii, jj, bsize[ib],
                                                      bsize[jb]);
                                    int ir = ROWMAJOR(jj, ii, bsize[jb],
                                                      bsize[ib]);
                                    double tmp = A_value_r[il];
                                    A_value_l[il] = A_value_r[ir];
                                    A_value_l[ir] = tmp;
                                }
                        }
                        else
                        {
                            for (int ii = 0; ii < bsize[ib]; ii++)
                                for (int jj = 0; jj < bsize[jb]; jj++)
                                {
                                    int il = ROWMAJOR(ii, jj, bsize[ib],
                                                      bsize[jb]);
                                    int ir = ROWMAJOR(jj, ii, bsize[jb],
                                                      bsize[ib]);
                                    double tmp = A_value_l[il];
                                    A_value_l[il] = A_value_r[ir];
                                    A_value_r[ir] = tmp;
                                }
                        }
                        exchangeDone = 1;
                        break;
                    }
                }
                assert(exchangeDone);
                // If no match add to end of row
//                if (!exchangeDone)
//                {
//                    int jind = A_nnzb[ind];
//                    {
//                        A_index[ROWMAJOR(ind, jind, N, M)] = i;
//                        A_value[ROWMAJOR(ind, jind, N, M)] =
//                            A_value[ROWMAJOR(i, j, N, M)];
//                        A_nnz[ind]++;
//                        A_nnz[i]--;
//                    }
//                }
            }
        }
    }

}
