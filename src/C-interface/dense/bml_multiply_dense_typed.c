#include "../../internal-blas/bml_gemm.h"
#include "../bml_logger.h"
#include "../../typed.h"
#include "bml_multiply.h"
#include "bml_trace.h"
#include "bml_parallel.h"
#include "bml_allocate.h"
#include "bml_multiply_dense.h"
#include "bml_trace_dense.h"
#include "bml_allocate_dense.h"
#include "bml_types.h"
#include "bml_types_dense.h"

#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/** Matrix multiply.
 *
 * \f$ C \leftarrow \alpha \, A \, B + \beta C \f$
 *
 * \ingroup multiply_group
 *
 * Note that the internal storage order of the dense matrix is
 * row-major while the standard implementation of gemm assumes column
 * major order. We therefore use a storage order transpose to the gemm
 * implementation. Instead of multiplying \f$ A \, B f$ we will
 * instead calculate \f$ B \, A \f$ and get:
 *
 * \f$ C^{T} \leftarrow \alpha \, B^{T} \, A^{T} + \beta C^{T} \f$
 *
 * \f$ C \f$ will then become:
 *
 * \f$ C \leftarrow \alpha \, \left( B^{T} \, A^{T} \right)^{T} + \beta C \f$
 * \f$ C \leftarrow \alpha \, A \, B + \beta C \f$
 *
 * \param A Matrix A
 * \param B Matrix B
 * \param C Matrix C
 * \param alpha Scalar factor multiplied by A * B
 * \param beta Scalar factor multiplied by C
 */
void TYPED_FUNC(
    bml_multiply_dense) (
    const bml_matrix_dense_t * A,
    const bml_matrix_dense_t * B,
    bml_matrix_dense_t * C,
    const double alpha,
    const double beta)
{
    REAL_T alpha_ = (REAL_T) alpha;
    REAL_T beta_ = (REAL_T) beta;

    TYPED_FUNC(bml_gemm) ("N", "N", &A->N, &A->N, &A->N, &alpha_, B->matrix,
                          &A->N, A->matrix, &A->N, &beta_, C->matrix, &A->N);
}

/** Matrix multiply.
 *
 * X2 = X * X
 *
 *  \ingroup multiply_group
 *
 *  \param X Matrix X
 *  \param X2 MatrixX2
 */
void *TYPED_FUNC(
    bml_multiply_x2_dense) (
    const bml_matrix_dense_t * X,
    bml_matrix_dense_t * X2)
{
    double *trace = bml_allocate_memory(sizeof(double) * 2);

    trace[0] = TYPED_FUNC(bml_trace_dense) (X);
    TYPED_FUNC(bml_multiply_dense) (X, X, X2, 1.0, 0.0);
    trace[1] = TYPED_FUNC(bml_trace_dense) (X2);

    return trace;
}

/** Matrix multiply.
 *
 * \f$ C \leftarrow A \, B \f$
 *
 * \ingroup multiply_group
 *
 * \param A Matrix A
 * \param B Matrix B
 * \param C Matrix C
 */
void TYPED_FUNC(
    bml_multiply_AB_dense) (
    const bml_matrix_dense_t * A,
    const bml_matrix_dense_t * B,
    bml_matrix_dense_t * C)
{
    TYPED_FUNC(bml_multiply_dense) (A, B, C, 1.0, 0.0);
}

/** Matrix multiply.
 *
 * This routine is provided for completeness.
 *
 * C = A * B
 *
 *  \ingroup multiply_group
 *
 *  \param A Matrix A
 *  \param B Matrix B
 *  \param C Matrix C
 */
void TYPED_FUNC(
    bml_multiply_adjust_AB_dense) (
    const bml_matrix_dense_t * A,
    const bml_matrix_dense_t * B,
    bml_matrix_dense_t * C)
{
    REAL_T alpha = (REAL_T) 1.0;
    REAL_T beta = (REAL_T) 0.0;

    TYPED_FUNC(bml_gemm) ("T", "T", &A->N, &A->N, &A->N, &alpha, A->matrix,
                          &A->N, B->matrix, &A->N, &beta, C->matrix, &A->N);
}
