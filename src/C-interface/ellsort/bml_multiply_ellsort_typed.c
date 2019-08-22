#include "../../macros.h"
#include "../../typed.h"
#include "../bml_add.h"
#include "../bml_allocate.h"
#include "../bml_logger.h"
#include "../bml_multiply.h"
#include "../bml_parallel.h"
#include "../bml_types.h"
#include "bml_add_ellsort.h"
#include "bml_allocate_ellsort.h"
#include "bml_multiply_ellsort.h"
#include "bml_types_ellsort.h"

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
    bml_multiply_ellsort) (
    bml_matrix_ellsort_t * A,
    bml_matrix_ellsort_t * B,
    bml_matrix_ellsort_t * C,
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
        TYPED_FUNC(bml_multiply_x2_ellsort) (A, C, threshold);
    }
    else
    {
        bml_matrix_dimension_t matrix_dimension = { A->N, A->N, A->M };
        bml_matrix_ellsort_t *A2 =
            TYPED_FUNC(bml_noinit_matrix_ellsort) (matrix_dimension,
                                                   A->distribution_mode);

        if (A != NULL && A == B)
        {
            TYPED_FUNC(bml_multiply_x2_ellsort) (A, A2, threshold);
        }
        else
        {
            TYPED_FUNC(bml_multiply_AB_ellsort) (A, B, A2, threshold);
        }

#ifdef DO_MPI
        if (bml_getNRanks() > 1 && A2->distribution_mode == distributed)
        {
            bml_allGatherVParallel(A2);
        }
#endif

        TYPED_FUNC(bml_add_ellsort) (C, A2, beta, alpha, threshold);

        bml_deallocate_ellsort(A2);
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
    bml_multiply_x2_ellsort) (
    bml_matrix_ellsort_t * X,
    bml_matrix_ellsort_t * X2,
    double threshold)
{
    int *X_localRowMin = X->domain->localRowMin;
    int *X_localRowMax = X->domain->localRowMax;

    int X_N = X->N;
    int X_M = X->M;
    int *X_index = X->index;
    int *X_nnz = X->nnz;

    int X2_N = X2->N;
    int X2_M = X2->M;
    int *X2_index = X2->index;
    int *X2_nnz = X2->nnz;

    REAL_T traceX = 0.0;
    REAL_T traceX2 = 0.0;
    REAL_T *X_value = (REAL_T *) X->value;
    REAL_T *X2_value = (REAL_T *) X2->value;

    double *trace = bml_allocate_memory(sizeof(double) * 2);

    int myRank = bml_getMyRank();

#if !(defined(__IBMC__) || defined(__ibmxl__))
    int ix[X_N], jx[X_N];
    REAL_T x[X_N];

    memset(ix, 0, X_N * sizeof(int));
    memset(jx, 0, X_N * sizeof(int));
    memset(x, 0.0, X_N * sizeof(REAL_T));
#endif

#if defined(__IBMC__) || defined(__ibmxl__)
#pragma omp parallel for                             \
  shared(X_N, X_M, X_index, X_nnz, X_value, myRank)  \
  shared(X2_N, X2_M, X2_index, X2_nnz, X2_value)     \
  shared(X_localRowMin, X_localRowMax)               \
  reduction(+: traceX, traceX2)
#else
#pragma omp parallel for                             \
  shared(X_N, X_M, X_index, X_nnz, X_value, myRank)  \
  shared(X2_N, X2_M, X2_index, X2_nnz, X2_value)     \
  shared(X_localRowMin, X_localRowMax)               \
  firstprivate(ix, jx, x)                            \
  reduction(+: traceX, traceX2)
#endif

    //for (int i = 0; i < X_N; i++)       // CALCULATES THRESHOLDED X^2
    for (int i = X_localRowMin[myRank]; i < X_localRowMax[myRank]; i++) // CALCULATES THRESHOLDED X^2
    {

#if defined(__IBMC__) || defined(__ibmxl__)
        int ix[X_N], jx[X_N];
        REAL_T x[X_N];

        memset(ix, 0, X_N * sizeof(int));
#endif

        int l = 0;
        for (int jp = 0; jp < X_nnz[i]; jp++)
        {
            REAL_T a = X_value[ROWMAJOR(i, jp, X_N, X_M)];
            int j = X_index[ROWMAJOR(i, jp, X_N, X_M)];
            if (j == i)
            {
                traceX = traceX + a;
            }
            for (int kp = 0; kp < X_nnz[j]; kp++)
            {
                int k = X_index[ROWMAJOR(j, kp, X_N, X_M)];
                if (ix[k] == 0)
                {
                    x[k] = 0.0;
                    //X2_index[ROWMAJOR(i, l, N, M)] = k;
                    jx[l] = k;
                    ix[k] = i + 1;
                    l++;
                }
                // TEMPORARY STORAGE VECTOR LENGTH FULL N
                x[k] = x[k] + a * X_value[ROWMAJOR(j, kp, X_N, X_M)];
            }
        }

        // Check for number of non-zeroes per row exceeded
        if (l > X2_M)
        {
            LOG_ERROR("Number of non-zeroes per row > M, Increase M\n");
        }

        int ll = 0;
        for (int j = 0; j < l; j++)
        {
            //int jp = X2_index[ROWMAJOR(i, j, N, M)];
            int jp = jx[j];
            REAL_T xtmp = x[jp];
            if (jp == i)
            {
                traceX2 = traceX2 + xtmp;
                X2_value[ROWMAJOR(i, ll, X2_N, X2_M)] = xtmp;
                X2_index[ROWMAJOR(i, ll, X2_N, X2_M)] = jp;
                ll++;
            }
            else if (is_above_threshold(xtmp, threshold))
            {
                X2_value[ROWMAJOR(i, ll, X2_N, X2_M)] = xtmp;
                X2_index[ROWMAJOR(i, ll, X2_N, X2_M)] = jp;
                ll++;
            }
            ix[jp] = 0;
            x[jp] = 0.0;
        }
        X2_nnz[i] = ll;
    }

    trace[0] = traceX;
    trace[1] = traceX2;

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
    bml_multiply_AB_ellsort) (
    bml_matrix_ellsort_t * A,
    bml_matrix_ellsort_t * B,
    bml_matrix_ellsort_t * C,
    double threshold)
{
    int A_N = A->N;
    int A_M = A->M;
    int *A_nnz = A->nnz;
    int *A_index = A->index;
    int *A_localRowMin = A->domain->localRowMin;
    int *A_localRowMax = A->domain->localRowMax;

    int B_N = B->N;
    int B_M = B->M;
    int *B_nnz = B->nnz;
    int *B_index = B->index;

    int C_N = C->N;
    int C_M = C->M;
    int *C_nnz = C->nnz;
    int *C_index = C->index;

    REAL_T *A_value = (REAL_T *) A->value;
    REAL_T *B_value = (REAL_T *) B->value;
    REAL_T *C_value = (REAL_T *) C->value;

    int myRank = bml_getMyRank();

#if !(defined(__IBMC__) || defined(__ibmxl__))
    int ix[C->N], jx[C->N];
    REAL_T x[C->N];

    memset(ix, 0, C->N * sizeof(int));
    memset(jx, 0, C->N * sizeof(int));
    memset(x, 0.0, C->N * sizeof(REAL_T));
#endif

#if defined(__IBMC__) || defined(__ibmxl__)
#pragma omp parallel for                     \
  shared(A_N, A_M, A_nnz, A_index, A_value)  \
  shared(A_localRowMin, A_localRowMax)       \
  shared(B_N, B_M, B_nnz, B_index, B_value)  \
  shared(C_N, C_M, C_nnz, C_index, C_value)  \
  shared(myRank)
#else
#pragma omp parallel for                     \
  shared(A_N, A_M, A_nnz, A_index, A_value)  \
  shared(A_localRowMin, A_localRowMax)       \
  shared(B_N, B_M, B_nnz, B_index, B_value)  \
  shared(C_N, C_M, C_nnz, C_index, C_value)  \
  shared(myRank)                             \
  firstprivate(ix, jx, x)
#endif

    //for (int i = 0; i < A_N; i++)
    for (int i = A_localRowMin[myRank]; i < A_localRowMax[myRank]; i++)
    {

#if defined(__IBMC__) || defined(__ibmxl__)
        int ix[C_N], jx[C_N];
        REAL_T x[C_N];

        memset(ix, 0, C_N * sizeof(int));
#endif

        int l = 0;
        for (int jp = 0; jp < A_nnz[i]; jp++)
        {
            REAL_T a = A_value[ROWMAJOR(i, jp, A_N, A_M)];
            int j = A_index[ROWMAJOR(i, jp, A_N, A_M)];

            for (int kp = 0; kp < B_nnz[j]; kp++)
            {
                int k = B_index[ROWMAJOR(j, kp, B_N, B_M)];
                if (ix[k] == 0)
                {
                    x[k] = 0.0;
                    //C_index[ROWMAJOR(i, l, N, M)] = k;
                    jx[l] = k;
                    ix[k] = i + 1;
                    l++;
                }
                // TEMPORARY STORAGE VECTOR LENGTH FULL N
                x[k] = x[k] + a * B_value[ROWMAJOR(j, kp, B_N, B_M)];
            }
        }

        // Check for number of non-zeroes per row exceeded
        if (l > C_M)
        {
            LOG_ERROR("Number of non-zeroes per row > M, Increase M\n");
        }

        int ll = 0;
        for (int j = 0; j < l; j++)
        {
            //int jp = C_index[ROWMAJOR(i, j, N, M)];
            int jp = jx[j];
            REAL_T xtmp = x[jp];
            if (jp == i)
            {
                C_value[ROWMAJOR(i, ll, C_N, C_M)] = xtmp;
                C_index[ROWMAJOR(i, ll, C_N, C_M)] = jp;
                ll++;
            }
            else if (is_above_threshold(xtmp, threshold))
            {
                C_value[ROWMAJOR(i, ll, C_N, C_M)] = xtmp;
                C_index[ROWMAJOR(i, ll, C_N, C_M)] = jp;
                ll++;
            }
            ix[jp] = 0;
            x[jp] = 0.0;
        }
        C_nnz[i] = ll;
    }
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
    bml_multiply_adjust_AB_ellsort) (
    bml_matrix_ellsort_t * A,
    bml_matrix_ellsort_t * B,
    bml_matrix_ellsort_t * C,
    double threshold)
{
    int A_N = A->N;
    int A_M = A->M;
    int *A_nnz = A->nnz;
    int *A_index = A->index;
    int *A_localRowMin = A->domain->localRowMin;
    int *A_localRowMax = A->domain->localRowMax;

    int B_N = B->N;
    int B_M = B->M;
    int *B_nnz = B->nnz;
    int *B_index = B->index;

    int C_N = C->N;
    int C_M = C->M;
    int *C_nnz = C->nnz;
    int *C_index = C->index;

    int aflag = 1;

    REAL_T *A_value = (REAL_T *) A->value;
    REAL_T *B_value = (REAL_T *) B->value;
    REAL_T *C_value = (REAL_T *) C->value;

    REAL_T adjust_threshold = (REAL_T) threshold;

    int myRank = bml_getMyRank();

#if !(defined(__IBMC__) || defined(__ibmxl__))
    int ix[C->N], jx[C->N];
    REAL_T x[C->N];

    memset(ix, 0, C->N * sizeof(int));
    memset(jx, 0, C->N * sizeof(int));
    memset(x, 0.0, C->N * sizeof(REAL_T));
#endif

    while (aflag > 0)
    {
        aflag = 0;

#if defined(__IBMC__) || defined(__ibmxl__)
#pragma omp parallel for                     \
  shared(A_N, A_M, A_nnz, A_index, A_value)  \
  shared(A_localRowMin, A_localRowMax)       \
  shared(B_N, B_M, B_nnz, B_index, B_value)  \
  shared(C_N, C_M, C_nnz, C_index, C_value)  \
  shared(adjust_threshold, myRank)           \
  reduction(+:aflag)
#else
#pragma omp parallel for                     \
  shared(A_N, A_M, A_nnz, A_index, A_value)  \
  shared(A_localRowMin, A_localRowMax)       \
  shared(B_N, B_M, B_nnz, B_index, B_value)  \
  shared(C_N, C_M, C_nnz, C_index, C_value)  \
  shared(adjust_threshold, myRank)           \
  firstprivate(ix, jx, x)                    \
  reduction(+:aflag)
#endif

        //for (int i = 0; i < A_N; i++)
        for (int i = A_localRowMin[myRank]; i < A_localRowMax[myRank]; i++)
        {

#if defined(__IBMC__) || defined(__ibmxl__)
            int ix[C_N], jx[C_N];
            REAL_T x[C_N];

            memset(ix, 0, C_N * sizeof(int));
#endif

            int l = 0;
            for (int jp = 0; jp < A_nnz[i]; jp++)
            {
                REAL_T a = A_value[ROWMAJOR(i, jp, A_N, A_M)];
                int j = A_index[ROWMAJOR(i, jp, A_N, A_M)];

                for (int kp = 0; kp < B_nnz[j]; kp++)
                {
                    int k = B_index[ROWMAJOR(j, kp, B_N, B_M)];
                    if (ix[k] == 0)
                    {
                        x[k] = 0.0;
                        jx[l] = k;
                        ix[k] = i + 1;
                        l++;
                    }
                    // TEMPORARY STORAGE VECTOR LENGTH FULL N
                    x[k] = x[k] + a * B_value[ROWMAJOR(j, kp, B_N, B_M)];
                }
            }

            // Check for number of non-zeroes per row exceeded
            // Need to adjust threshold
            if (l > C_M)
            {
                aflag = 1;
            }

            int ll = 0;
            for (int j = 0; j < l; j++)
            {
                //int jp = C_index[ROWMAJOR(i, j, N, M)];
                int jp = jx[j];
                REAL_T xtmp = x[jp];
                // Diagonal elements are saved in first column
                if (jp == i)
                {
                    C_value[ROWMAJOR(i, ll, C_N, C_M)] = xtmp;
                    C_index[ROWMAJOR(i, ll, C_N, C_M)] = jp;
                    ll++;
                }
                else if (is_above_threshold(xtmp, adjust_threshold))
                {
                    C_value[ROWMAJOR(i, ll, C_N, C_M)] = xtmp;
                    C_index[ROWMAJOR(i, ll, C_N, C_M)] = jp;
                    ll++;
                }
                ix[jp] = 0;
                x[jp] = 0.0;
            }
            C_nnz[i] = ll;
        }

        adjust_threshold *= (REAL_T) 2.0;
    }
}
