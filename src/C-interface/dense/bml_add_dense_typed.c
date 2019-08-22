/* Needs to be included before #include <complex.h>. */
#ifdef BML_USE_MAGMA
#include "magma_v2.h"
#endif

#include "../../macros.h"
#include "../../typed.h"
#include "../blas.h"
#include "../bml_add.h"
#include "../bml_allocate.h"
#include "../bml_logger.h"
#include "../bml_parallel.h"
#include "../bml_scale.h"
#include "../bml_types.h"
#include "bml_add_dense.h"
#include "bml_allocate_dense.h"
#include "bml_copy_dense.h"
#include "bml_scale_dense.h"
#include "bml_types_dense.h"

#include <complex.h>
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
 */
void TYPED_FUNC(
    bml_add_dense) (
    bml_matrix_dense_t * A,
    bml_matrix_dense_t * B,
    double alpha,
    double beta)
{
    int myRank = bml_getMyRank();

    int nElems = B->domain->localRowExtent[myRank] * B->N;
    int startIndex = B->domain->localDispl[myRank];
    int inc = 1;

#ifdef BML_USE_MAGMA
    nElems = B->N * B->ld;
    MAGMA_T alpha_ = MAGMACOMPLEX(MAKE) (alpha, 0.);
    MAGMA_T beta_ = MAGMACOMPLEX(MAKE) (beta, 0.);
    MAGMA(scal) (nElems, alpha_, A->matrix, inc, A->queue);
    MAGMA(axpy) (nElems, beta_, B->matrix, inc,
                 A->matrix + startIndex, inc, A->queue);
#else
    REAL_T alpha_ = alpha;
    REAL_T beta_ = beta;

#ifdef NOBLAS
    LOG_ERROR("No BLAS library");
#else
    C_BLAS(SCAL) (&nElems, &alpha_, A->matrix + startIndex, &inc);
    C_BLAS(AXPY) (&nElems, &beta_, B->matrix + startIndex, &inc,
                  A->matrix + startIndex, &inc);
#endif

#endif
}

/** Matrix addition and calculate TrNorm.
 *
 * \f$ A = \alpha A + \beta B \f$
 *
 * \ingroup add_group
 *
 * \param A Matrix A
 * \param B Matrix B
 * \param alpha Scalar factor multiplied by A
 * \param beta Scalar factor multiplied by B
 */
double TYPED_FUNC(
    bml_add_norm_dense) (
    bml_matrix_dense_t * A,
    bml_matrix_dense_t * B,
    double alpha,
    double beta)
{
    double trnorm = 0.0;
    REAL_T *B_matrix = (REAL_T *) B->matrix;
    int myRank = bml_getMyRank();
    int N = A->N;

    int *A_localRowMin = A->domain->localRowMin;
    int *A_localRowMax = A->domain->localRowMax;

#pragma omp parallel for                                \
  shared(B_matrix, A_localRowMin, A_localRowMax)        \
  shared(N, myRank)                                     \
  reduction(+:trnorm)
    //for (int i = 0; i < N * N; i++)
    for (int i = A_localRowMin[myRank] * N; i < A_localRowMax[myRank] * N;
         i++)
    {
        trnorm += B_matrix[i] * B_matrix[i];
    }

    TYPED_FUNC(bml_add_dense) (A, B, alpha, beta);

    return trnorm;
}

/** Matrix addition.
 *
 * \f$ A = A + \beta \mathrm{Id} \f$
 *
 *  \ingroup add_group
 *
 *  \param A Matrix A
 *  \param beta Scalar factor multiplied by I
 */
void TYPED_FUNC(
    bml_add_identity_dense) (
    bml_matrix_dense_t * A,
    double beta)
{
    int N = A->N;
#if BML_USE_MAGMA
    MAGMA_T *A_matrix = (MAGMA_T *) A->matrix;
    MAGMA_T beta_ = MAGMACOMPLEX(MAKE) (beta, 0.);
    bml_matrix_dense_t *B =
        TYPED_FUNC(bml_identity_matrix_dense) (N, sequential);
    MAGMABLAS(geadd) (N, N, beta_, (MAGMA_T *) B->matrix, B->ld,
                      A_matrix, A->ld, A->queue);
    bml_deallocate_dense(B);
#else
    REAL_T *A_matrix = (REAL_T *) A->matrix;
    REAL_T beta_ = beta;
    int *A_localRowMin = A->domain->localRowMin;
    int *A_localRowMax = A->domain->localRowMax;
    int myRank = bml_getMyRank();

#pragma omp parallel for                                \
  shared(A_matrix, A_localRowMin, A_localRowMax)        \
  shared(N, myRank, beta_)
    //for (int i = 0; i < N; i++)
    for (int i = A_localRowMin[myRank]; i < A_localRowMax[myRank]; i++)
    {
        A_matrix[ROWMAJOR(i, i, N, N)] += beta_;
    }
#endif
}

/** Matrix addition.
 *
 * \f$ A = alpha A + \beta \mathrm{Id} \f$
 *
 *  \ingroup add_group
 *
 *  \param A Matrix A
 *  \param alpha Scalar factor multiplied by A
 *  \param beta Scalar factor multiplied by I
 */
void TYPED_FUNC(
    bml_scale_add_identity_dense) (
    bml_matrix_dense_t * A,
    double alpha,
    double beta)
{
    REAL_T _alpha = (REAL_T) alpha;

    bml_scale_inplace_dense(&_alpha, A);
    bml_add_identity_dense(A, beta);
}
