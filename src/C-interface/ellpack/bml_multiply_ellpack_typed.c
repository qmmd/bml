#include "../../macros.h"
#include "../../typed.h"
#include "../bml_add.h"
#include "../bml_allocate.h"
#include "../bml_logger.h"
#include "../bml_multiply.h"
#include "../bml_parallel.h"
#include "../bml_types.h"
#include "bml_add_ellpack.h"
#include "bml_allocate_ellpack.h"
#include "bml_multiply_ellpack.h"
#include "bml_types_ellpack.h"

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

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
    bml_multiply_ellpack) (
    bml_matrix_ellpack_t * A,
    bml_matrix_ellpack_t * B,
    bml_matrix_ellpack_t * C,
    double alpha,
    double beta,
    double threshold)
{
    double ONE = 1.0;
    double ZERO = 0.0;

    void *trace = NULL;

    if (A == NULL || B == NULL)
    {
        LOG_ERROR("Either matrix A or B are NULL\n");
    }

    if (A == B && alpha == ONE && beta == ZERO)
    {
        trace = TYPED_FUNC(bml_multiply_x2_ellpack) (A, C, threshold);
    }
    else
    {
        bml_matrix_dimension_t matrix_dimension = { C->N, C->N, C->M };
        bml_matrix_ellpack_t *A2 =
            TYPED_FUNC(bml_noinit_matrix_ellpack) (matrix_dimension,
                                                   A->distribution_mode);

        if (A != NULL && A == B)
        {
            trace = TYPED_FUNC(bml_multiply_x2_ellpack) (A, A2, threshold);
        }
        else
        {
            TYPED_FUNC(bml_multiply_AB_ellpack) (A, B, A2, threshold);
        }

#ifdef DO_MPI
        if (bml_getNRanks() > 1 && A2->distribution_mode == distributed)
        {
            bml_allGatherVParallel(A2);
        }
#endif

        TYPED_FUNC(bml_add_ellpack) (C, A2, beta, alpha, threshold);

        bml_deallocate_ellpack(A2);
    }
    bml_free_memory(trace);
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
    bml_multiply_x2_ellpack) (
    bml_matrix_ellpack_t * X,
    bml_matrix_ellpack_t * X2,
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
    int rowMin = X_localRowMin[myRank];
    int rowMax = X_localRowMax[myRank];

#if !(defined(__IBMC__) || defined(__ibmxl__) || (defined(USE_OMP_OFFLOAD) && (defined(INTEL_SDK) || defined(CRAY_SDK))))
    int ix[X_N], jx[X_N];
    REAL_T x[X_N];

    memset(ix, 0, X_N * sizeof(int));
    memset(jx, 0, X_N * sizeof(int));
    memset(x, 0.0, X_N * sizeof(REAL_T));
#endif

#if defined(USE_OMP_OFFLOAD) && (defined(INTEL_SDK) || defined(CRAY_SDK))
    int num_chunks = MIN(OFFLOAD_NUM_CHUNKS, rowMax - rowMin + 1);

    int all_ix[X_N * num_chunks], all_jx[X_N * num_chunks];
    REAL_T all_x[X_N * num_chunks];

    memset(all_ix, 0, X_N * num_chunks * sizeof(int));
    memset(all_jx, 0, X_N * num_chunks * sizeof(int));
    memset(all_x, 0.0, X_N * num_chunks * sizeof(REAL_T));

#pragma omp target map(to:all_ix[0:X_N*num_chunks],all_jx[0:X_N*num_chunks],all_x[0:X_N*num_chunks])

#endif

#if defined (USE_OMP_OFFLOAD)
#if defined(INTEL_SDK) || defined(CRAY_SDK)
#pragma omp teams distribute parallel for	\
    shared(X_N, X_M, X_index, X_nnz, X_value)  \
    shared(X2_N, X2_M, X2_index, X2_nnz, X2_value)     \
    shared(rowMin, rowMax)                             \
    reduction(+: traceX, traceX2)
    for (int chunk = 0; chunk < num_chunks; chunk++)
    {
        int *ix, *jx;
        REAL_T *x;

        ix = &all_ix[chunk * X_N];
        jx = &all_jx[chunk * X_N];
        x = &all_x[chunk * X_N];

#else
#pragma omp target teams distribute parallel for                               \
    shared(X_N, X_M, X_index, X_nnz, X_value)  \
    shared(X2_N, X2_M, X2_index, X2_nnz, X2_value)     \
    shared(rowMin, rowMax)                             \
    firstprivate(ix,jx, x)                             \
    reduction(+: traceX, traceX2)
#endif
#else
#if defined(__IBMC__) || defined(__ibmxl__)
#pragma omp parallel for                               \
    shared(X_N, X_M, X_index, X_nnz, X_value)  \
    shared(X2_N, X2_M, X2_index, X2_nnz, X2_value)     \
    shared(rowMin, rowMax)                             \
    reduction(+: traceX, traceX2)
#else
#if !(defined(CRAY_SDK) || defined(INTEL_SDK))
#pragma vector aligned
#endif
#pragma omp parallel for                               \
    shared(X_N, X_M, X_index, X_nnz, X_value)  \
    shared(X2_N, X2_M, X2_index, X2_nnz, X2_value)     \
    shared(rowMin, rowMax)                             \
    firstprivate(ix,jx, x)                             \
    reduction(+: traceX, traceX2)
#endif
#endif

#if defined(USE_OMP_OFFLOAD) && (defined(INTEL_SDK) || defined(CRAY_SDK))
    for (int i = rowMin + chunk; i < rowMax; i = i + num_chunks)
#else
    for (int i = rowMin; i < rowMax; i++)
#endif
    {

#if defined(__IBMC__) || defined(__ibmxl__)
        int ix[X_N], jx[X_N];
        REAL_T x[X_N];

        memset(ix, 0, X_N * sizeof(int));
#endif

#ifdef INTEL_OPT
        __assume_aligned(X_nnz, MALLOC_ALIGNMENT);
        __assume_aligned(X_index, MALLOC_ALIGNMENT);
        __assume_aligned(X_value, MALLOC_ALIGNMENT);
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
#ifndef USE_OMP_OFFLOAD
            LOG_ERROR("Number of non-zeroes per row > M, Increase M\n");
#endif
        }

#ifdef INTEL_OPT
        __assume_aligned(X2_nnz, MALLOC_ALIGNMENT);
        __assume_aligned(X2_index, MALLOC_ALIGNMENT);
        __assume_aligned(X2_value, MALLOC_ALIGNMENT);
#endif
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

#if defined(USE_OMP_OFFLOAD) && (defined(INTEL_SDK) || defined(CRAY_SDK))
}
#endif

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
    bml_multiply_AB_ellpack) (
    bml_matrix_ellpack_t * A,
    bml_matrix_ellpack_t * B,
    bml_matrix_ellpack_t * C,
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
    int rowMin = A_localRowMin[myRank];
    int rowMax = A_localRowMax[myRank];

#if !(defined(__IBMC__) || defined(__ibmxl__) || (defined(USE_OMP_OFFLOAD) && (defined(INTEL_SDK) || defined(CRAY_SDK))))
    int ix[C->N], jx[C->N];
    REAL_T x[C->N];

    memset(ix, 0, C->N * sizeof(int));
    memset(jx, 0, C->N * sizeof(int));
    memset(x, 0.0, C->N * sizeof(REAL_T));
#endif

#if defined(USE_OMP_OFFLOAD) && (defined(INTEL_SDK) || defined(CRAY_SDK))
    int num_chunks = MIN(OFFLOAD_NUM_CHUNKS, rowMax - rowMin + 1);

    int all_ix[C_N * num_chunks], all_jx[C_N * num_chunks];
    REAL_T all_x[C_N * num_chunks];

    memset(all_ix, 0, C_N * num_chunks * sizeof(int));
    memset(all_jx, 0, C_N * num_chunks * sizeof(int));
    memset(all_x, 0.0, C_N * num_chunks * sizeof(REAL_T));

#pragma omp target map(to:all_ix[0:C_N*num_chunks],all_jx[0:C_N*num_chunks],all_x[0:C_N*num_chunks])

#endif

#if defined (USE_OMP_OFFLOAD)
#if defined(INTEL_SDK) || defined(CRAY_SDK)
#pragma omp teams distribute parallel for \
    shared(A_N, A_M, A_nnz, A_index, A_value)  \
    shared(A_localRowMin, A_localRowMax)       \
    shared(B_N, B_M, B_nnz, B_index, B_value)  \
    shared(C_N, C_M, C_nnz, C_index, C_value)
    for (int chunk = 0; chunk < num_chunks; chunk++)
    {
        int *ix, *jx;
        REAL_T *x;

        ix = &all_ix[chunk * C_N];
        jx = &all_jx[chunk * C_N];
        x = &all_x[chunk * C_N];

#else
#pragma omp target teams distribute parallel for \
    shared(A_N, A_M, A_nnz, A_index, A_value)  \
    shared(A_localRowMin, A_localRowMax)       \
    shared(B_N, B_M, B_nnz, B_index, B_value)  \
    shared(C_N, C_M, C_nnz, C_index, C_value)  \
    firstprivate(ix, jx, x)
#endif
#else
#if defined(__IBMC__) || defined(__ibmxl__)
#pragma omp parallel for                       \
    shared(A_N, A_M, A_nnz, A_index, A_value)  \
    shared(A_localRowMin, A_localRowMax)       \
    shared(B_N, B_M, B_nnz, B_index, B_value)  \
    shared(C_N, C_M, C_nnz, C_index, C_value)
#else
#pragma omp parallel for                       \
    shared(A_N, A_M, A_nnz, A_index, A_value)  \
    shared(A_localRowMin, A_localRowMax)       \
    shared(B_N, B_M, B_nnz, B_index, B_value)  \
    shared(C_N, C_M, C_nnz, C_index, C_value)  \
    firstprivate(ix, jx, x)
#endif
#endif
    //for (int i = 0; i < A_N; i++)
#if defined(USE_OMP_OFFLOAD) && (defined(INTEL_SDK) || defined(CRAY_SDK))
    for (int i = rowMin + chunk; i < rowMax; i = i + num_chunks)
#else
    for (int i = rowMin; i < rowMax; i++)
#endif
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
#ifndef USE_OMP_OFFLOAD
            LOG_ERROR("Number of non-zeroes per row > M, Increase M\n");
#endif
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
#if defined(USE_OMP_OFFLOAD) && (defined(INTEL_SDK) || defined(CRAY_SDK))
}
#endif
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
    bml_multiply_adjust_AB_ellpack) (
    bml_matrix_ellpack_t * A,
    bml_matrix_ellpack_t * B,
    bml_matrix_ellpack_t * C,
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
    int rowMin = A_localRowMin[myRank];
    int rowMax = A_localRowMax[myRank];

#ifdef USE_OMP_OFFLOAD
#pragma omp target update from(A_nnz[:A_N], A_index[:A_N*A_M], A_value[:A_N*A_M])
#pragma omp target update from(B_nnz[:B_N], B_index[:B_N*B_M], B_value[:B_N*B_M])
#endif

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
#pragma omp parallel for                       \
    shared(A_N, A_M, A_nnz, A_index, A_value)  \
    shared(B_N, B_M, B_nnz, B_index, B_value)  \
    shared(C_N, C_M, C_nnz, C_index, C_value)  \
    shared(adjust_threshold)           \
    reduction(+:aflag)
#else
#pragma omp parallel for                       \
    shared(A_N, A_M, A_nnz, A_index, A_value)  \
    shared(B_N, B_M, B_nnz, B_index, B_value)  \
    shared(C_N, C_M, C_nnz, C_index, C_value)  \
    shared(adjust_threshold)           \
    firstprivate(ix, jx, x)                    \
    reduction(+:aflag)
#endif

        //for (int i = 0; i < A_N; i++)
        for (int i = rowMin; i < rowMax; i++)
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
#ifdef USE_OMP_OFFLOAD
#pragma omp target update to(C_nnz[:C_N], C_index[:C_N*C_M], C_value[:C_N*C_M])
#endif
}
