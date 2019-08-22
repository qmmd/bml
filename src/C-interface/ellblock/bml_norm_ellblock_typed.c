#include "../../macros.h"
#include "../../typed.h"
#include "../bml_norm.h"
#include "../bml_parallel.h"
#include "../bml_types.h"
#include "bml_norm_ellblock.h"
#include "bml_types_ellblock.h"
#include "bml_utilities_ellblock.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/** Calculate the sum of squares of the elements of a matrix.
 *
 *  \ingroup norm_group
 *
 *  \param A The matrix A
 *  \return The sum of squares of A
 */
double TYPED_FUNC(
    bml_sum_squares_ellblock) (
    bml_matrix_ellblock_t * A)
{
    int NB = A->NB;
    int MB = A->MB;

    int *A_nnzb = A->nnzb;
    int *A_indexb = A->indexb;
    int *bsize = A->bsize;

    REAL_T sum = 0.0;
    REAL_T **A_ptr_value = (REAL_T **) A->ptr_value;

    for (int ib = 0; ib < NB; ib++)
    {
        for (int jp = 0; jp < A_nnzb[ib]; jp++)
        {
            int ind = ROWMAJOR(ib, jp, NB, MB);
            int jb = A_indexb[ind];
            REAL_T *xval = A_ptr_value[ind];
            for (int ii = 0; ii < bsize[ib]; ii++)
                for (int jj = 0; jj < bsize[jb]; jj++)
                {
                    int index = ROWMAJOR(ii, jj, bsize[ib], bsize[jb]);
                    sum += xval[index] * xval[index];
                }
        }
    }

    return (double) REAL_PART(sum);
}

/** Calculate the sum of squares of the elements of \alpha A + \beta B.
 *
 *  \ingroup norm_group
 *
 *  \param A The matrix A
 *  \param B The matrix B
 *  \param alpha Multiplier for A
 *  \param beta Multiplier for B
 *  \pram threshold Threshold
 *  \return The sum of squares of \alpha A + \beta B
 */
double TYPED_FUNC(
    bml_sum_squares2_ellblock) (
    bml_matrix_ellblock_t * A,
    bml_matrix_ellblock_t * B,
    double alpha,
    double beta,
    double threshold)
{
    int NB = A->NB;
    int MB = A->MB;

    int *A_indexb = A->indexb;
    int *A_nnzb = A->nnzb;
    int *bsize = A->bsize;
    int *B_indexb = B->indexb;
    int *B_nnzb = B->nnzb;

    REAL_T sum = 0.0;
    REAL_T **A_ptr_value = (REAL_T **) A->ptr_value;
    REAL_T **B_ptr_value = (REAL_T **) B->ptr_value;

    REAL_T alpha_ = (REAL_T) alpha;
    REAL_T beta_ = (REAL_T) beta;

    int maxbsize = 0;
    for (int ib = 0; ib < NB; ib++)
        maxbsize = MAX(maxbsize, bsize[ib]);
    int maxbsize2 = maxbsize * maxbsize;
    REAL_T *y_ptr[NB];
    for (int ib = 0; ib < NB; ib++)
        y_ptr[ib] = calloc(maxbsize2, sizeof(REAL_T));

    int ix[NB], jjb[NB];

    memset(ix, 0, NB * sizeof(int));
    memset(jjb, 0, NB * sizeof(int));

    for (int ib = 0; ib < NB; ib++)
    {
        int lb = 0;
        for (int jp = 0; jp < A_nnzb[ib]; jp++)
        {
            int ind = ROWMAJOR(ib, jp, NB, MB);
            int jb = A_indexb[ind];
            int nelements = bsize[ib] * bsize[jb];
            if (ix[jb] == 0)
            {
                memset(y_ptr[jb], 0, nelements * sizeof(REAL_T));
                ix[jb] = ib + 1;
                jjb[lb] = jb;
                lb++;
            }
            REAL_T *y_value = y_ptr[jb];
            REAL_T *A_value = A_ptr_value[ind];
            for (int ii = 0; ii < bsize[ib]; ii++)
                for (int jj = 0; jj < bsize[jb]; jj++)
                {
                    int index = ROWMAJOR(ii, jj, bsize[ib], bsize[jb]);
                    y_value[index] += alpha_ * A_value[index];
                }
        }

        for (int jp = 0; jp < B_nnzb[ib]; jp++)
        {
            int ind = ROWMAJOR(ib, jp, NB, MB);
            int jb = B_indexb[ind];
            int nelements = bsize[ib] * bsize[jb];
            if (ix[jb] == 0)
            {
                memset(y_ptr[jb], 0, nelements * sizeof(REAL_T));
                ix[jb] = ib + 1;
                jjb[lb] = jb;
                lb++;
            }
            REAL_T *y_value = y_ptr[jb];
            REAL_T *B_value = B_ptr_value[ind];
            for (int ii = 0; ii < bsize[ib]; ii++)
                for (int jj = 0; jj < bsize[jb]; jj++)
                {
                    int index = ROWMAJOR(ii, jj, bsize[ib], bsize[jb]);
                    y_value[index] += beta_ * B_value[index];
                }
        }

        for (int jp = 0; jp < lb; jp++)
        {
            double normx = TYPED_FUNC(bml_sum_squares)
                (y_ptr[jjb[jp]], bsize[ib], bsize[jp], bsize[jp]);

            if (normx > threshold * threshold)
                sum += normx;

            ix[jjb[jp]] = 0;
            memset(y_ptr[jjb[jp]], 0, maxbsize2 * sizeof(REAL_T));
            jjb[jp] = 0;
        }
    }

    for (int ib = 0; ib < NB; ib++)
        free(y_ptr[ib]);
    return (double) REAL_PART(sum);
}

/** Calculate the Frobenius norm of matrix A.
 *
 *  \ingroup norm_group
 *
 *  \param A The matrix A
 *  \return The Frobenius norm of A
 */
double TYPED_FUNC(
    bml_fnorm_ellblock) (
    bml_matrix_ellblock_t * A)
{
    double fnorm = TYPED_FUNC(bml_sum_squares_ellblock) (A);
#ifdef DO_MPI
    if (bml_getNRanks() > 1 && A->distribution_mode == distributed)
    {
        bml_sumRealReduce(&fnorm);
    }
#endif
    fnorm = sqrt(fnorm);

    return (double) REAL_PART(fnorm);
}

/** Calculate the Frobenius norm of 2 matrices.
 *
 *  \ingroup norm_group
 *
 *  \param A The matrix A
 *  \param B The matrix B
 *  \return The Frobenius norm of A-B
 */
double TYPED_FUNC(
    bml_fnorm2_ellblock) (
    bml_matrix_ellblock_t * A,
    bml_matrix_ellblock_t * B)
{
    int NB = A->NB;
    int MB = A->MB;
    double fnorm = 0.0;
    REAL_T *rvalue;

    int *A_nnzb = (int *) A->nnzb;
    int *A_indexb = (int *) A->indexb;
    int *bsize = A->bsize;

    REAL_T **A_ptr_value = (REAL_T **) A->ptr_value;
    int *B_nnzb = (int *) B->nnzb;
    int *B_indexb = (int *) B->indexb;
    REAL_T **B_ptr_value = (REAL_T **) B->ptr_value;

    REAL_T temp;

    int maxbsize = 0;
    for (int ib = 0; ib < NB; ib++)
        maxbsize = MAX(maxbsize, bsize[ib]);
    int maxbsize2 = maxbsize * maxbsize;
    REAL_T *zero_block = calloc(maxbsize2, sizeof(REAL_T));

    for (int ib = 0; ib < NB; ib++)
    {
        for (int jp = 0; jp < A_nnzb[ib]; jp++)
        {
            int ind = ROWMAJOR(ib, jp, NB, MB);
            for (int kp = 0; kp < B_nnzb[ib]; kp++)
            {
                if (A_indexb[ROWMAJOR(ib, jp, NB, MB)] ==
                    B_indexb[ROWMAJOR(ib, kp, NB, MB)])
                {
                    rvalue = B_ptr_value[ROWMAJOR(ib, kp, NB, MB)];
                    break;
                }
                rvalue = zero_block;
            }

            REAL_T *A_value = A_ptr_value[ind];
            int jb = A_indexb[ind];
            for (int ii = 0; ii < bsize[ib]; ii++)
                for (int jj = 0; jj < bsize[jb]; jj++)
                {
                    int index = ROWMAJOR(ii, jj, bsize[ib], bsize[jb]);
                    temp = A_value[index] - rvalue[index];
                    fnorm += temp * temp;
                }
        }

        for (int jp = 0; jp < B_nnzb[ib]; jp++)
        {
            for (int kp = 0; kp < A_nnzb[ib]; kp++)
            {
                if (A_indexb[ROWMAJOR(ib, kp, NB, MB)] ==
                    B_indexb[ROWMAJOR(ib, jp, NB, MB)])
                {
                    rvalue = A_ptr_value[ROWMAJOR(ib, kp, NB, MB)];
                    break;
                }
                rvalue = NULL;
            }

            if (rvalue == NULL)
            {
                int ind = ROWMAJOR(ib, jp, NB, MB);
                int jb = A_indexb[ind];
                REAL_T *B_value = B_ptr_value[ind];
                for (int ii = 0; ii < bsize[ib]; ii++)
                    for (int jj = 0; jj < bsize[jb]; jj++)
                    {
                        int index = ROWMAJOR(ii, jj, bsize[ib], bsize[jb]);
                        temp = B_value[index];
                        fnorm += temp * temp;
                    }
            }
        }
    }

#ifdef DO_MPI
    if (bml_getNRanks() > 1 && A->distribution_mode == distributed)
    {
        bml_sumRealReduce(&fnorm);
    }
#endif

    fnorm = sqrt(fnorm);

    free(zero_block);

    return (double) REAL_PART(fnorm);
}
