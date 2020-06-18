#include "../../macros.h"
#include "../../typed.h"
#include "../bml_add.h"
#include "../bml_allocate.h"
#include "../bml_logger.h"
#include "../bml_multiply.h"
#include "../bml_parallel.h"
#include "../bml_types.h"
#include "bml_add_csr.h"
#include "bml_allocate_csr.h"
#include "bml_multiply_csr.h"
#include "bml_types_csr.h"
#include "bml_setters_csr.h"

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
    bml_multiply_csr) (
    bml_matrix_csr_t * A,
    bml_matrix_csr_t * B,
    bml_matrix_csr_t * C,
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
        trace = TYPED_FUNC(bml_multiply_x2_csr) (A, C, threshold);
    }
    else
    {
        bml_matrix_dimension_t matrix_dimension = { C->N_, C->N_, C->NZMAX_ };
        bml_matrix_csr_t *A2 =
            TYPED_FUNC(bml_noinit_matrix_csr) (matrix_dimension,
                                               A->distribution_mode);

        if (A != NULL && A == B)
        {
            trace = TYPED_FUNC(bml_multiply_x2_csr) (A, A2, threshold);
        }
        else
        {
            TYPED_FUNC(bml_multiply_AB_csr) (A, B, A2, threshold);
        }

#ifdef DO_MPI
        if (bml_getNRanks() > 1 && A2->distribution_mode == distributed)
        {
            bml_allGatherVParallel(A2);
        }
#endif

        TYPED_FUNC(bml_add_csr) (C, A2, beta, alpha, threshold);

        bml_deallocate_csr(A2);
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
    bml_multiply_x2_csr) (
    bml_matrix_csr_t * X,
    bml_matrix_csr_t * X2,
    double threshold)
{
    int X_N = X->N_;

    REAL_T traceX = 0.0;
    REAL_T traceX2 = 0.0;

    double *trace = bml_allocate_memory(sizeof(double) * 2);

#if !(defined(__IBMC__) || defined(__ibmxl__))
    int ix[X_N], jx[X_N];
    REAL_T x[X_N];

    memset(ix, 0, X_N * sizeof(int));
    memset(jx, 0, X_N * sizeof(int));
    memset(x, 0.0, X_N * sizeof(REAL_T));
#endif

#if defined(__IBMC__) || defined(__ibmxl__)
#pragma omp parallel for                               \
    shared(X_N)  \
    reduction(+: traceX, traceX2)
#else
#if !(defined(CRAY_SDK) || defined(INTEL_SDK))
#pragma vector aligned
#endif
#pragma omp parallel for                               \
    shared(X_N)  \
    firstprivate(ix,jx, x)                             \
    reduction(+: traceX, traceX2)
#endif

    for (int i = 0; i < X_N; i++)       // CALCULATES THRESHOLDED X^2
    {

#if defined(__IBMC__) || defined(__ibmxl__)
        int ix[X_N], jx[X_N];
        REAL_T x[X_N];

        memset(ix, 0, X_N * sizeof(int));
#endif

        int *icols = X->data_[i]->cols_;
        REAL_T *ivals = (REAL_T *) X->data_[i]->vals_;
        const int innz = X->data_[i]->NNZ_;

        int l = 0;
        for (int ipos = 0; ipos < innz; ipos++)
        {
            REAL_T a = ivals[ipos];
            const int j = icols[ipos];
            if (j == i)
            {
                traceX = traceX + a;
            }
            const int jnnz = X->data_[j]->NNZ_;
            REAL_T *jvals = (REAL_T *) X->data_[j]->vals_;
            int *jcols = X->data_[j]->cols_;
            for (int jpos = 0; jpos < jnnz; jpos++)
            {
                const int k = jcols[jpos];
                if (ix[k] == 0)
                {
                    x[k] = 0.0;
                    jx[l] = k;
                    ix[k] = i + 1;
                    l++;
                }
                // TEMPORARY STORAGE VECTOR LENGTH FULL N
                x[k] = x[k] + a * jvals[jpos];
            }
        }

        // clear row
        TYPED_FUNC(csr_clear_row) (X2->data_[i]);
        for (int j = 0; j < l; j++)
        {
            int jp = jx[j];
            REAL_T xtmp = x[jp];
            if (jp == i)
            {
                traceX2 = traceX2 + xtmp;
                TYPED_FUNC(csr_set_row_element_new) (X2->data_[i], jp, &xtmp);
            }
            else if (is_above_threshold(xtmp, threshold))
            {
                TYPED_FUNC(csr_set_row_element_new) (X2->data_[i], jp, &xtmp);
            }
            // reset
            ix[jp] = 0;
            x[jp] = 0.0;
        }
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
    bml_multiply_AB_csr) (
    bml_matrix_csr_t * A,
    bml_matrix_csr_t * B,
    bml_matrix_csr_t * C,
    double threshold)
{
    const int A_N = A->N_;
    const int C_N = C->N_;

#if !(defined(__IBMC__) || defined(__ibmxl__))
    int ix[C_N], jx[C_N];
    REAL_T x[C_N];

    memset(ix, 0, C_N * sizeof(int));
    memset(jx, 0, C_N * sizeof(int));
    memset(x, 0.0, C_N * sizeof(REAL_T));
#endif

#if defined(__IBMC__) || defined(__ibmxl__)
#pragma omp parallel for                       \
    shared(A_N, C_N)
#else
#pragma omp parallel for                       \
    firstprivate(ix, jx, x)
#endif

    for (int i = 0; i < A_N; i++)
    {
#if defined(__IBMC__) || defined(__ibmxl__)
        int ix[C_N], jx[C_N];
        REAL_T x[C_N];

        memset(ix, 0, C_N * sizeof(int));
#endif
        int *acols = A->data_[i]->cols_;
        REAL_T *avals = (REAL_T *) A->data_[i]->vals_;
        const int annz = A->data_[i]->NNZ_;
        int l = 0;
        for (int pos = 0; pos < annz; pos++)
        {
            REAL_T a = avals[pos];
            const int j = acols[pos];

            const int bnnz = B->data_[j]->NNZ_;
            REAL_T *bvals = (REAL_T *) B->data_[j]->vals_;
            int *bcols = B->data_[j]->cols_;
            for (int bpos = 0; bpos < bnnz; bpos++)
            {
                const int k = bcols[bpos];
                if (ix[k] == 0)
                {
                    x[k] = 0.0;
                    jx[l] = k;
                    ix[k] = i + 1;
                    l++;
                }
                // TEMPORARY STORAGE VECTOR LENGTH FULL N
                x[k] = x[k] + a * bvals[bpos];
            }
        }
        // clear row
        TYPED_FUNC(csr_clear_row) (C->data_[i]);
        for (int j = 0; j < l; j++)
        {
            int jp = jx[j];
            REAL_T xtmp = x[jp];
            if (jp == i)
            {
                TYPED_FUNC(csr_set_row_element_new) (C->data_[i], jp, &xtmp);
            }
            else if (is_above_threshold(xtmp, threshold))
            {
                TYPED_FUNC(csr_set_row_element_new) (C->data_[i], jp, &xtmp);
            }
            // reset
            ix[jp] = 0;
            x[jp] = 0.0;
        }
    }
}

/* Not sure why we have this routine (for ellpack, ellblock, etc). Using default here,
* with no threshold adjustment.
*/
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
    bml_multiply_adjust_AB_csr) (
    bml_matrix_csr_t * A,
    bml_matrix_csr_t * B,
    bml_matrix_csr_t * C,
    double threshold)
{
    TYPED_FUNC(bml_multiply_AB_csr) (A, B, C, threshold);
}
