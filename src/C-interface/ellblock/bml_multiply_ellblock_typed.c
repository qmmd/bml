#include "../../internal-blas/bml_gemm.h"
#include "../../macros.h"
#include "../../typed.h"
#include "../bml_add.h"
#include "../bml_allocate.h"
#include "../bml_logger.h"
#include "../bml_multiply.h"
#include "../bml_multiply.h"
#include "../bml_parallel.h"
#include "../bml_types.h"
#include "bml_add_ellblock.h"
#include "bml_allocate_ellblock.h"
#include "bml_multiply_ellblock.h"
#include "bml_types_ellblock.h"
#include "bml_utilities_ellblock.h"

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/** Matrix multiply.
 *
 * \f$ C \leftarrow \alpha A \, B + \beta C \f$
 *
 * \ingroup multiply_group
 *
 * \param A Matrix A
 * \param B Matrix B
 * \param C Matrix C
 * \param alpha Scalar factor multiplied by A * B
 * \param beta Scalar factor multiplied by C
 * \param threshold Used for sparse multiply
 */
void TYPED_FUNC(
    bml_multiply_ellblock) (
    bml_matrix_ellblock_t * A,
    bml_matrix_ellblock_t * B,
    bml_matrix_ellblock_t * C,
    double alpha,
    double beta,
    double threshold)
{
    double ONE = 1.0;
    double ZERO = 0.0;

    if (A == NULL || B == NULL)
    {
        LOG_ERROR("Either matrix A or B are NULL\n");
    }

    if (A == B && alpha == ONE && beta == ZERO)
    {
        TYPED_FUNC(bml_multiply_x2_ellblock) (A, C, threshold);
    }
    else
    {
        bml_matrix_ellblock_t *A2 =
            TYPED_FUNC(bml_block_matrix_ellblock) (C->NB, C->MB, C->M,
                                                   C->bsize,
                                                   C->distribution_mode);

        if (A != NULL && A == B)
        {
            TYPED_FUNC(bml_multiply_x2_ellblock) (A, A2, threshold);
        }
        else
        {
            TYPED_FUNC(bml_multiply_AB_ellblock) (A, B, A2, threshold);
        }

#ifdef DO_MPI
        if (bml_getNRanks() > 1 && A2->distribution_mode == distributed)
        {
            bml_allGatherVParallel(A2);
        }
#endif

        TYPED_FUNC(bml_add_ellblock) (C, A2, beta, alpha, threshold);

        bml_deallocate_ellblock(A2);
    }
}

/** Matrix multiply.
 *
 * \f$ X^{2} \leftarrow X \, X \f$
 *
 * \ingroup multiply_group
 *
 * \param X Matrix X
 * \param X2 Matrix X2
 * \param threshold Used for sparse multiply
 */
void *TYPED_FUNC(
    bml_multiply_x2_ellblock) (
    bml_matrix_ellblock_t * X,
    bml_matrix_ellblock_t * X2,
    double threshold)
{
    int NB = X->NB;
    int MB = X->MB;
    int *X_indexb = X->indexb;
    int *X_nnzb = X->nnzb;
    int *bsize = X->bsize;

    int *X2_indexb = X2->indexb;
    int *X2_nnzb = X2->nnzb;

    int ix[NB], jx[NB];
    REAL_T *x_ptr[NB];

    REAL_T traceX = 0.0;
    REAL_T traceX2 = 0.0;
    REAL_T **X_ptr_value = (REAL_T **) X->ptr_value;
    REAL_T **X2_ptr_value = (REAL_T **) X2->ptr_value;

    double *trace = bml_allocate_memory(sizeof(double) * 2);

    memset(ix, 0, NB * sizeof(int));
    memset(jx, 0, NB * sizeof(int));

    int maxbsize = 0;
    for (int ib = 0; ib < NB; ib++)
        maxbsize = MAX(maxbsize, bsize[ib]);
    int maxbsize2 = maxbsize * maxbsize;
#ifdef _OPENMP
    const int nthreads = omp_get_max_threads();
#else
    const int nthreads = 1;
#endif
    REAL_T *x_ptr_storage = calloc(maxbsize2 * NB * nthreads, sizeof(REAL_T));

    char xptrset = 0;

#pragma omp parallel for                           \
    firstprivate(ix,jx, x_ptr, xptrset)            \
    reduction(+: traceX, traceX2)

    //loop over row blocks
    for (int ib = 0; ib < NB; ib++)
    {
        int lb = 0;
        if (!xptrset)
        {
#ifdef _OPENMP
            int offset = omp_get_thread_num() * maxbsize2 * NB;
#else
            int offset = 0;
#endif
            for (int i = 0; i < NB; i++)
                x_ptr[i] = &x_ptr_storage[offset + i * maxbsize2];
            xptrset = 1;
        }

        // loop over non-zero blocks in row block ib
        for (int jp = 0; jp < X_nnzb[ib]; jp++)
        {
            int ind = ROWMAJOR(ib, jp, NB, MB);
            REAL_T *X_value_left = X_ptr_value[ind];
            int jb = X_indexb[ind];
            if (jb == ib)
            {
                for (int kk = 0; kk < bsize[ib]; kk++)
                    traceX += X_value_left[kk * (1 + bsize[ib])];
            }
            for (int kp = 0; kp < X_nnzb[jb]; kp++)
            {
                int indk = ROWMAJOR(jb, kp, NB, MB);
                int kb = X_indexb[indk];
                if (ix[kb] == 0)
                {
                    memset(x_ptr[kb], 0.0, maxbsize2 * sizeof(REAL_T));
                    jx[lb] = kb;
                    ix[kb] = ib + 1;
                    lb++;
                }
                REAL_T *x = x_ptr[kb];
                REAL_T *X_value_right = X_ptr_value[indk];

                // multiply block ib,jb by block jb,kb
#ifdef NOBLAS
                for (int ii = 0; ii < bsize[ib]; ii++)
                    for (int jj = 0; jj < bsize[kb]; jj++)
                    {
                        int k = ROWMAJOR(ii, jj, bsize[ib], bsize[kb]);
                        for (int kk = 0; kk < bsize[jb]; kk++)
                        {
                            x[k] = x[k]
                                +
                                X_value_left[ROWMAJOR
                                             (ii, kk, bsize[ib],
                                              bsize[jb])] *
                                X_value_right[ROWMAJOR
                                              (kk, jj, bsize[jb], bsize[kb])];
                        }
                    }
#else
                REAL_T alpha = (REAL_T) 1.;
                REAL_T beta = (REAL_T) 1.;
                TYPED_FUNC(bml_xsmm_gemm) ("N", "N", &bsize[kb], &bsize[ib],
                                           &bsize[jb], &alpha, X_value_right,
                                           &bsize[kb], X_value_left,
                                           &bsize[jb], &beta, x, &bsize[kb]);
#endif
            }
        }

        // Check for number of non-zero blocks per block row exceeded
        if (lb > MB)
        {
            LOG_ERROR
                ("Number of non-zeroes blocks per row > MB, Increase MB\n");
        }

        int ll = 0;
        for (int jb = 0; jb < lb; jb++)
        {
            assert(ll < MB);
            int jp = jx[jb];
            REAL_T *xtmp = x_ptr[jp];
            double normx = TYPED_FUNC(bml_norm_inf_fast)
                (xtmp, bsize[ib] * bsize[jp]);
            if (jp == ib || (normx > threshold))
            {
                if (jp == ib)
                {
                    for (int kk = 0; kk < bsize[ib]; kk++)
                        traceX2 = traceX2 + xtmp[kk * (1 + bsize[jb])];
                }

                int nelements = bsize[ib] * bsize[jp];
                int ind = ROWMAJOR(ib, ll, NB, MB);
                assert(ind < NB * MB);
                if (X2_ptr_value[ind] == NULL)
                    X2_ptr_value[ind] =
                        TYPED_FUNC(bml_allocate_block_ellblock) (X2, ib,
                                                                 nelements);
                REAL_T *X2_value = X2_ptr_value[ind];
                assert(X2_value != NULL);
                memcpy(X2_value, xtmp, nelements * sizeof(REAL_T));
                X2_indexb[ind] = jp;
                ll++;
            }
            ix[jp] = 0;
            memset(x_ptr[jp], 0, maxbsize2 * sizeof(REAL_T));
        }
        X2_nnzb[ib] = ll;
    }

    trace[0] = traceX;
    trace[1] = traceX2;

    free(x_ptr_storage);

    return trace;
}

/** Matrix multiply.
 *
 * \f$ C \leftarrow B \, A \f$
 *
 * \ingroup multiply_group
 *
 * \param A Matrix A
 * \param B Matrix B
 * \param C Matrix C
 * \param threshold Used for sparse multiply
 */
void TYPED_FUNC(
    bml_multiply_AB_ellblock) (
    bml_matrix_ellblock_t * A,
    bml_matrix_ellblock_t * B,
    bml_matrix_ellblock_t * C,
    double threshold)
{
    assert(A->NB == B->NB);
    assert(A->NB == C->NB);

    int NB = A->NB;

    int *A_nnzb = A->nnzb;
    int *A_indexb = A->indexb;
    int *bsize = A->bsize;

    int *B_nnzb = B->nnzb;
    int *B_indexb = B->indexb;

    int *C_nnzb = C->nnzb;
    int *C_indexb = C->indexb;

    int ix[NB], jx[NB];
    REAL_T *x_ptr[NB];

    REAL_T **A_ptr_value = (REAL_T **) A->ptr_value;
    REAL_T **B_ptr_value = (REAL_T **) B->ptr_value;
    REAL_T **C_ptr_value = (REAL_T **) C->ptr_value;

    memset(ix, 0, NB * sizeof(int));
    memset(jx, 0, NB * sizeof(int));

    int maxbsize = 0;
    for (int ib = 0; ib < NB; ib++)
        maxbsize = MAX(maxbsize, bsize[ib]);
    int maxbsize2 = maxbsize * maxbsize;
#ifdef _OPENMP
    int nthreads = omp_get_max_threads();
#else
    int nthreads = 1;
#endif
    REAL_T *x_ptr_storage = calloc(maxbsize2 * NB * nthreads, sizeof(REAL_T));

    char xptrset = 0;

    //loop over row blocks
#pragma omp parallel for                       \
    firstprivate(ix, jx, x_ptr, xptrset)

    for (int ib = 0; ib < NB; ib++)
    {
        int lb = 0;
        if (!xptrset)
        {
#ifdef _OPENMP
            int offset = omp_get_thread_num() * maxbsize2 * NB;
#else
            int offset = 0;
#endif
            for (int i = 0; i < NB; i++)
                x_ptr[i] = &x_ptr_storage[offset + i * maxbsize2];
            xptrset = 1;
        }

        //loop over blocks in this block row "ib"
        for (int jp = 0; jp < A_nnzb[ib]; jp++)
        {
            int ind = ROWMAJOR(ib, jp, NB, A->MB);
            REAL_T *A_value = A_ptr_value[ind];
            int jb = A_indexb[ind];

            for (int kp = 0; kp < B_nnzb[jb]; kp++)
            {
                int kb = B_indexb[ROWMAJOR(jb, kp, NB, B->MB)];
                //compute column block "kb" of result
                if (ix[kb] == 0)
                {
                    memset(x_ptr[kb], 0, maxbsize2 * sizeof(REAL_T));
                    jx[lb] = kb;
                    ix[kb] = ib + 1;
                    lb++;
                }
                REAL_T *x = x_ptr[kb];
                REAL_T *B_value = B_ptr_value[ROWMAJOR(jb, kp, NB, B->MB)];
                for (int ii = 0; ii < bsize[ib]; ii++)
                    for (int jj = 0; jj < bsize[kb]; jj++)
                    {
                        int k = ROWMAJOR(ii, jj, bsize[ib], bsize[kb]);
                        for (int kk = 0; kk < bsize[jb]; kk++)
                        {
                            x[k] = x[k]
                                +
                                A_value[ROWMAJOR
                                        (ii, kk, bsize[ib],
                                         bsize[jb])] * B_value[ROWMAJOR(kk,
                                                                        jj,
                                                                        bsize
                                                                        [jb],
                                                                        bsize
                                                                        [kb])];
                        }
                    }
            }
        }

        // Check for number of non-zeroes per row exceeded
        if (lb > C->MB)
        {
            LOG_ERROR("Number of non-zeroes per row > M, Increase M\n");
        }

        int ll = 0;
        for (int jb = 0; jb < lb; jb++)
        {
            int jp = jx[jb];
            REAL_T *xtmp = x_ptr[jp];
            double normx = TYPED_FUNC(bml_norm_inf)
                (xtmp, bsize[ib], bsize[jp], bsize[jp]);

            if (jp == ib || (normx > threshold))
            {
                int nelements = bsize[ib] * bsize[jp];
                int ind = ROWMAJOR(ib, ll, NB, C->MB);
                C_ptr_value[ind]
                    =
                    TYPED_FUNC(bml_allocate_block_ellblock) (C, ib,
                                                             nelements);
                memcpy(C_ptr_value[ind], xtmp, nelements * sizeof(REAL_T));
                C_indexb[ind] = jp;
                ll++;
            }
            ix[jp] = 0;
            memset(x_ptr[jp], 0, maxbsize2 * sizeof(REAL_T));
        }
        C_nnzb[ib] = ll;
    }

    free(x_ptr_storage);
}

/** Matrix multiply with threshold adjustment.
 *
 * \f$ C \leftarrow B \, A \f$
 *
 * \ingroup multiply_group
 *
 * \param A Matrix A
 * \param B Matrix B
 * \param C Matrix C
 * \param threshold Used for sparse multiply
 */
void TYPED_FUNC(
    bml_multiply_adjust_AB_ellblock) (
    bml_matrix_ellblock_t * A,
    bml_matrix_ellblock_t * B,
    bml_matrix_ellblock_t * C,
    double threshold)
{
    int NB = A->NB;
    int MB = A->MB;
    int *A_nnzb = A->nnzb;
    int *A_indexb = A->indexb;
    int *bsize = A->bsize;

    int *B_nnzb = B->nnzb;
    int *B_indexb = B->indexb;

    int *C_nnzb = C->nnzb;
    int *C_indexb = C->indexb;

    int ix[NB], jx[NB];
    int aflag = 1;
    REAL_T *x_ptr[NB];

    REAL_T **A_ptr_value = (REAL_T **) A->ptr_value;
    REAL_T **B_ptr_value = (REAL_T **) B->ptr_value;
    REAL_T **C_ptr_value = (REAL_T **) C->ptr_value;

    double adjust_threshold = threshold;

    memset(ix, 0, NB * sizeof(int));
    memset(jx, 0, NB * sizeof(int));

    int maxbsize = 0;
    for (int ib = 0; ib < NB; ib++)
        maxbsize = MAX(maxbsize, bsize[ib]);
    int maxbsize2 = maxbsize * maxbsize;
    for (int ib = 0; ib < NB; ib++)
        x_ptr[ib] = calloc(maxbsize2, sizeof(REAL_T));

    while (aflag > 0)
    {
        aflag = 0;

        for (int ib = 0; ib < NB; ib++)
        {
            int l = 0;
            for (int jp = 0; jp < A_nnzb[ib]; jp++)
            {
                REAL_T *a = A_ptr_value[ROWMAJOR(ib, jp, NB, MB)];
                int jb = A_indexb[ROWMAJOR(ib, jp, NB, MB)];

                for (int kp = 0; kp < B_nnzb[jb]; kp++)
                {
                    int kb = B_indexb[ROWMAJOR(jb, kp, NB, MB)];
                    if (ix[kb] == 0)
                    {
                        memset(x_ptr[kb], 0, maxbsize2);
                        jx[l] = kb;
                        ix[kb] = ib + 1;
                        l++;
                    }
                    // TEMPORARY STORAGE VECTOR LENGTH FULL N
                    REAL_T *b = B_ptr_value[ROWMAJOR(jb, kp, NB, MB)];
                    REAL_T *x = x_ptr[kb];
                    for (int ii = 0; ii < bsize[jb]; ii++)
                        for (int jj = 0; jj < bsize[kb]; jj++)
                        {
                            int k = ii * bsize[kb] + jj;
                            x[k] =
                                x[k] +
                                a[ROWMAJOR(ii, jj, bsize[ib], bsize[jb])] *
                                b[ROWMAJOR(jj, ii, bsize[jb], bsize[kb])];
                        }
                }
            }

            // Check for number of non-zeroes per row exceeded
            // Need to adjust threshold
            if (l > MB)
            {
                aflag = 1;
            }

            int ll = 0;
            for (int jb = 0; jb < l; jb++)
            {
                int jp = jx[jb];
                REAL_T *xtmp = x_ptr[jp];
                double normx = TYPED_FUNC(bml_norm_inf)
                    (xtmp, bsize[ib], bsize[jp], bsize[jp]);
                // Diagonal elements are saved in first column
                if (jp == ib)
                {
                    memcpy(C_ptr_value[ROWMAJOR(ib, ll, NB, MB)], xtmp,
                           bsize[ib] * bsize[ll] * sizeof(REAL_T));
                    C_indexb[ROWMAJOR(ib, ll, NB, MB)] = jp;
                    ll++;
                }
                else if (normx > adjust_threshold)
                {
                    memcpy(C_ptr_value[ROWMAJOR(ib, ll, NB, MB)], xtmp,
                           bsize[ib] * bsize[ll] * sizeof(REAL_T));
                    C_indexb[ROWMAJOR(ib, ll, NB, MB)] = jp;
                    ll++;
                }
                ix[jp] = 0;
                memset(x_ptr[jp], 0, maxbsize2 * sizeof(REAL_T));
            }
            C_nnzb[ib] = ll;
        }

        adjust_threshold *= 2.0;
    }
}
