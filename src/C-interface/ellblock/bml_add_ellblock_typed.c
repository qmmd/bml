#include "../../macros.h"
#include "../../typed.h"
#include "../bml_add.h"
#include "../bml_allocate.h"
#include "../bml_parallel.h"
#include "../bml_types.h"
#include "bml_add_ellblock.h"
#include "bml_allocate_ellblock.h"
#include "bml_types_ellblock.h"
#include "bml_utilities_ellblock.h"

#include <assert.h>
#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/** Matrix addition.
 *
 * \f$ A = \alpha A + \beta B \f$
 *
 * \ingroup add_group
 *
 * \param A Matrix A
 * \param B Matrix B
 * \param alpha Scalar factor multiplied by A
 * \param beta Scalar factor multiplied by B
 * \param threshold Threshold for matrix addition
 */
void TYPED_FUNC(
    bml_add_ellblock) (
    bml_matrix_ellblock_t * A,
    bml_matrix_ellblock_t * B,
    double alpha,
    double beta,
    double threshold)
{
    assert(A->NB == B->NB);
    assert(A->bsize[0] == B->bsize[0]);

    int NB = A->NB;
    int MB = A->MB;
    int ix[NB], jx[NB];

    int *A_nnzb = A->nnzb;
    int *A_indexb = A->indexb;

    int *B_nnzb = B->nnzb;
    int *B_indexb = B->indexb;

    int *bsize = A->bsize;

    REAL_T *x_ptr[NB];
    REAL_T **A_ptr_value = (REAL_T **) A->ptr_value;
    REAL_T **B_ptr_value = (REAL_T **) B->ptr_value;

    memset(ix, 0, NB * sizeof(int));
    memset(jx, 0, NB * sizeof(int));

    int maxbsize = 0;
    for (int ib = 0; ib < NB; ib++)
        maxbsize = MAX(maxbsize, bsize[ib]);
    int maxbsize2 = maxbsize * maxbsize;
    for (int ib = 0; ib < NB; ib++)
        x_ptr[ib] = calloc(maxbsize2, sizeof(REAL_T));

    for (int ib = 0; ib < NB; ib++)
    {
        int l = 0;
        if (alpha > (double) 0.0 || alpha < (double) 0.0)
            for (int jp = 0; jp < A_nnzb[ib]; jp++)
            {
                int ind = ROWMAJOR(ib, jp, NB, MB);
                int jb = A_indexb[ind];
                REAL_T *x = x_ptr[jb];
                int nelements = bsize[ib] * bsize[jb];
                if (ix[jb] == 0)
                {
                    for (int kk = 0; kk < nelements; kk++)
                    {
                        x[kk] = 0.0;
                    }
                    ix[jb] = ib + 1;
                    jx[l] = jb;
                    l++;
                }
                REAL_T *A_value = A_ptr_value[ind];
                for (int kk = 0; kk < nelements; kk++)
                {
                    x[kk] = x[kk] + alpha * A_value[kk];
                }
            }

        if (beta > (double) 0.0 || beta < (double) 0.0)
            for (int jp = 0; jp < B_nnzb[ib]; jp++)
            {
                int ind = ROWMAJOR(ib, jp, NB, MB);
                int jb = B_indexb[ind];
                REAL_T *x = x_ptr[jb];
                int nelements = bsize[ib] * bsize[jb];
                if (ix[jb] == 0)
                {
                    for (int kk = 0; kk < nelements; kk++)
                    {
                        x[kk] = 0.0;
                    }
                    ix[jb] = ib + 1;
                    jx[l] = jb;
                    l++;
                }
                REAL_T *B_value = B_ptr_value[ind];
                for (int kk = 0; kk < nelements; kk++)
                {
                    x[kk] = x[kk] + beta * B_value[kk];
                }
            }
        A_nnzb[ib] = l;

        int ll = 0;
        for (int jp = 0; jp < l; jp++)
        {
            int jb = jx[jp];
            REAL_T *x = x_ptr[jb];
            REAL_T normx = TYPED_FUNC(bml_norm_inf)
                (x, bsize[ib], bsize[jb], bsize[jb]);
            int nelements = bsize[ib] * bsize[jb];
            if (is_above_threshold(normx, threshold))
            {
                int kb = ROWMAJOR(ib, ll, NB, MB);
                if (A_ptr_value[kb] == NULL)
                {
                    A_ptr_value[kb]
                        =
                        bml_noinit_allocate_memory(nelements *
                                                   sizeof(REAL_T));
                }
                for (int kk = 0; kk < nelements; kk++)
                {
                    A_ptr_value[kb][kk] = x[kk];
                }
                A_indexb[kb] = jb;
                ll++;
            }
            memset(x, 0.0, nelements * sizeof(REAL_T));
            ix[jb] = 0;
        }
        A_nnzb[ib] = ll;
    }

    for (int ib = 0; ib < NB; ib++)
        free(x_ptr[ib]);
}

/** Matrix addition.
 *
 * \f$ A = \alpha A + \beta B \f$
 *
 * \ingroup add_group
 *
 * \param A Matrix A
 * \param B Matrix B
 * \param alpha Scalar factor multiplied by A
 * \param beta Scalar factor multiplied by B
 * \param threshold Threshold for matrix addition
 */
double TYPED_FUNC(
    bml_add_norm_ellblock) (
    bml_matrix_ellblock_t * A,
    bml_matrix_ellblock_t * B,
    double alpha,
    double beta,
    double threshold)
{
    int NB = A->NB;
    int MB = A->MB;

    int ix[NB], jx[NB];

    int *A_nnzb = A->nnzb;
    int *A_indexb = A->indexb;

    int *B_nnzb = B->nnzb;
    int *B_indexb = B->indexb;

    int *bsize = A->bsize;

    REAL_T *x[NB];
    REAL_T **A_ptr_value = (REAL_T **) A->ptr_value;
    REAL_T **B_ptr_value = (REAL_T **) B->ptr_value;

    memset(ix, 0, NB * sizeof(int));
    memset(jx, 0, NB * sizeof(int));
    x[0] = (REAL_T *) calloc(BMAXSIZE * NB, sizeof(REAL_T));
    for (int i = 1; i < NB; i++)
    {
        x[i] = x[0] + i * BMAXSIZE;
    }

    REAL_T *y[NB];
    y[0] = (REAL_T *) calloc(BMAXSIZE * NB, sizeof(REAL_T));
    for (int i = 1; i < NB; i++)
    {
        y[i] = y[0] + i * BMAXSIZE;
    }


    double trnorm = 0.0;

    for (int ib = 0; ib < NB; ib++)
    {
        int l = 0;
        for (int jp = 0; jp < A_nnzb[ib]; jp++)
        {
            int ind = ROWMAJOR(ib, jp, NB, MB);
            int jb = A_indexb[ind];
            if (ix[jb] == 0)
            {
                for (int ii = 0; ii < bsize[ib]; ii++)
                    for (int jj = 0; jj < bsize[jb]; jj++)
                    {
                        int kk = ROWMAJOR(ii, jj, bsize[ib], bsize[jb]);
                        x[jb][kk] = 0.0;
                        y[jb][kk] = 0.0;
                    }
                ix[jb] = ib + 1;
                jx[l] = jb;
                l++;
            }
            for (int ii = 0; ii < bsize[ib]; ii++)
                for (int jj = 0; jj < bsize[jb]; jj++)
                {
                    int kk = ROWMAJOR(ii, jj, bsize[ib], bsize[jb]);
                    x[jb][kk] = x[jb][kk] + alpha * A_ptr_value[ind][kk];
                    y[jb][kk] = y[jb][kk] + A_ptr_value[ind][kk];
                }
        }

        for (int jp = 0; jp < B_nnzb[ib]; jp++)
        {
            int ind = ROWMAJOR(ib, jp, NB, MB);
            int jb = B_indexb[ind];
            if (ix[jb] == 0)
            {
                for (int ii = 0; ii < bsize[ib]; ii++)
                    for (int jj = 0; jj < bsize[jb]; jj++)
                    {
                        int kk = ROWMAJOR(ii, jj, bsize[ib], bsize[jb]);
                        x[jb][kk] = 0.0;
                        y[jb][kk] = 0.0;
                    }
                ix[jb] = ib + 1;
                jx[l] = jb;
                l++;
            }
            REAL_T *B_value = B_ptr_value[ind];
            for (int ii = 0; ii < bsize[ib]; ii++)
                for (int jj = 0; jj < bsize[jb]; jj++)
                {
                    int kk = ROWMAJOR(ii, jj, bsize[ib], bsize[jb]);
                    x[jb][kk] = x[jb][kk] + beta * B_ptr_value[ind][kk];
                    y[jb][kk] = y[jb][kk] - B_value[kk];
                }
        }
        A_nnzb[ib] = l;

        int ll = 0;
        for (int jp = 0; jp < l; jp++)
        {
            int jb = jx[jp];
            REAL_T *xTmp = x[jb];
            REAL_T normx = TYPED_FUNC(bml_norm_inf)
                (x[jb], bsize[ib], bsize[jb], bsize[jb]);
            REAL_T normy = TYPED_FUNC(bml_sum_squares)
                (y[jb], bsize[ib], bsize[jb], bsize[jb]);
            trnorm += normy * normy;
            if (is_above_threshold(normx, threshold))
            {
                int kb = ROWMAJOR(ib, ll, NB, MB);
                for (int ii = 0; ii < bsize[ib]; ii++)
                    for (int jj = 0; jj < bsize[jb]; jj++)
                    {
                        int kk = ROWMAJOR(ii, jj, bsize[ib], bsize[jb]);
                        A_ptr_value[kb][kk] = xTmp[kk];
                    }

                A_indexb[kb] = jb;
                ll++;
            }
            memset(x[jb], 0.0, bsize[ib] * bsize[jb] * sizeof(REAL_T));
            memset(y[jb], 0.0, bsize[ib] * bsize[jb] * sizeof(REAL_T));
            ix[jb] = 0;
        }
        A_nnzb[ib] = ll;
    }

    return trnorm;
}

/** Matrix addition.
 *
 *  \f$ A \leftarrow A + beta * I \f$
 *
 *  \ingroup add_group
 *
 *  \param A Matrix A
 *  \param beta Scalar factor multiplied by I
 *  \param threshold Threshold for matrix addition
 */
void TYPED_FUNC(
    bml_add_identity_ellblock) (
    bml_matrix_ellblock_t * A,
    double beta,
    double threshold)
{
    REAL_T alpha = (REAL_T) 1.0;

    bml_matrix_ellblock_t *Id =
        TYPED_FUNC(bml_identity_matrix_ellblock) (A->N, A->M,
                                                  A->distribution_mode);

    TYPED_FUNC(bml_add_ellblock) (A, Id, alpha, beta, threshold);

    bml_deallocate_ellblock(Id);
}

/** Matrix addition.
 *
 *  A = alpha * A + beta * I
 *
 *  \ingroup add_group
 *
 *  \param A Matrix A
 *  \param alpha Scalar factor multiplied by A
 *  \param beta Scalar factor multiplied by I
 *  \param threshold Threshold for matrix addition
 */
void TYPED_FUNC(
    bml_scale_add_identity_ellblock) (
    bml_matrix_ellblock_t * A,
    double alpha,
    double beta,
    double threshold)
{
    bml_matrix_ellblock_t *Id =
        TYPED_FUNC(bml_identity_matrix_ellblock) (A->N, A->M,
                                                  A->distribution_mode);

    TYPED_FUNC(bml_add_ellblock) (A, Id, alpha, beta, threshold);

    bml_deallocate_ellblock(Id);
}
