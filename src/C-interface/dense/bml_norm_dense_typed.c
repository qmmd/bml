#include "../../macros.h"
#include "../blas.h"
#include "../../typed.h"
#include "bml_norm.h"
#include "bml_norm_dense.h"
#include "bml_types.h"
#include "bml_types_dense.h"
#include "bml_parallel.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/** Calculate the sum of squares of all the elements of a matrix.
 *
 *  \ingroup norm_group
 *
 *  \param A The matrix
 *  \return The sum of squares of A
 */
double TYPED_FUNC(
    bml_sum_squares_dense) (
    const bml_matrix_dense_t * A)
{
    int N = A->N;

    REAL_T sum = 0.0;
    REAL_T *A_matrix = A->matrix;

    int *A_localRowMin = A->domain->localRowMin;
    int *A_localRowMax = A->domain->localRowMax;

    int myRank = bml_getMyRank();

#pragma omp parallel for default(none) \
    shared(N, A_matrix) \
    shared(A_localRowMin, A_localRowMax, myRank) \
    reduction(+:sum)
    for (int i = A_localRowMin[myRank] * N; i < A_localRowMax[myRank] * N;
         i++)
    {
        sum += A_matrix[i] * A_matrix[i];
    }

    return (double) REAL_PART(sum);
}

/** Calculate the sum of squares of all the core elements of a submatrix.
 *
 *  \ingroup norm_group
 *
 *  \param A The matrix
 *  \param core_size Number of core rows
 *  \return The sum of squares of A
 */
double TYPED_FUNC(
    bml_sum_squares_submatrix_dense) (
    const bml_matrix_dense_t * A,
    const int core_size)
{
    int N = A->N;

    REAL_T sum = 0.0;
    REAL_T *A_matrix = A->matrix;

#pragma omp parallel for default(none) \
    shared(N, A_matrix) \
    reduction(+:sum)
    for (int i = 0; i < core_size * N; i++)
    {
        sum += A_matrix[i] * A_matrix[i];
    }

/*
    for (int i = 0; i < core_size; i++)
    {
        //for (int j = 0; j < core_size; j++)
        for (int j = 0; j < N; j++)
        {
            REAL_T value = A_matrix[ROWMAJOR(i, j, N, N)];
            sum += value * value;
        }
    }
*/

    return (double) REAL_PART(sum);
}

/** Calculate the sum of squares of all elements of \alpha A + \beta B.
 *
 *  \ingroup norm_group
 *
 *  \param A The matrix A
 *  \param B The matrix B
 *  \param alpha Multiplier for A
 *  \param beta Multiplier for B
 *  \param threshold Threshold
 *  \return The sum of squares of all elements of \alpha A + \beta BB
 */
double TYPED_FUNC(
    bml_sum_squares2_dense) (
    const bml_matrix_dense_t * A,
    const bml_matrix_dense_t * B,
    const double alpha,
    const double beta,
    const double threshold)
{
    int N = A->N;

    REAL_T sum = 0.0;
    REAL_T *A_matrix = A->matrix;
    REAL_T *B_matrix = B->matrix;

    int *A_localRowMin = A->domain->localRowMin;
    int *A_localRowMax = A->domain->localRowMax;

    REAL_T alpha_ = (REAL_T) alpha;
    REAL_T beta_ = (REAL_T) beta;

    int myRank = bml_getMyRank();

#pragma omp parallel for \
    default(none) \
    shared(alpha_, beta_) \
    shared(N, A_matrix, B_matrix) \
    shared(A_localRowMin, A_localRowMax, myRank) \
    reduction(+:sum)
    //for (int i = 0; i < N * N; i++)
    for (int i = A_localRowMin[myRank] * N; i < A_localRowMax[myRank] * N;
         i++)
    {
        REAL_T temp = alpha_ * A_matrix[i] + beta_ * B_matrix[i];
        if (ABS(temp) > threshold)
            sum += temp * temp;
    }

    return (double) REAL_PART(sum);
}

/** Calculate the Frobenius norm of a matrix.
 *
 *  \ingroup norm_group
 *
 *  \param A The matrix
 *  \return The Frobenius norm of A
 */
double TYPED_FUNC(
    bml_fnorm_dense) (
    const bml_matrix_dense_t * A)
{
    double sum = 0.0;
    sum = TYPED_FUNC(bml_sum_squares_dense) (A);

#ifdef DO_MPI
    if (bml_getNRanks() > 1 && A->distribution_mode == distributed)
    {
        bml_sumRealReduce(&sum);
    }
#endif

    return (double) REAL_PART(sqrt(sum));
}

/** Calculate the Frobenius norm of 2 matrices.
 *
 *  \ingroup norm_group
 *
 *  \param A The matrix A
 *  \param A The matrix B
 *  \return The Frobenius norm of A-B
 */
double TYPED_FUNC(
    bml_fnorm2_dense) (
    const bml_matrix_dense_t * A,
    const bml_matrix_dense_t * B)
{
    int N = A->N;

    REAL_T *A_matrix = (REAL_T *) A->matrix;
    REAL_T *B_matrix = (REAL_T *) B->matrix;

    int *A_localRowMin = A->domain->localRowMin;
    int *A_localRowMax = A->domain->localRowMax;

    double fnorm = 0.0;
    REAL_T temp;

    int myRank = bml_getMyRank();

#pragma omp parallel for \
    default(none) \
    shared(temp) \
    shared(N, A_matrix, B_matrix) \
    shared(A_localRowMin, A_localRowMax, myRank) \
    reduction(+:fnorm)
    //for (int i = 0; i < N*N;i++)
    for (int i = A_localRowMin[myRank] * N; i < A_localRowMax[myRank] * N;
         i++)
    {
        temp = A_matrix[i] - B_matrix[i];
        fnorm += temp * temp;
    }

#ifdef DO_MPI
    if (bml_getNRanks() > 1 && A->distribution_mode == distributed)
    {
        bml_sumRealReduce(&fnorm);
    }
#endif

    return (double) REAL_PART(sqrt(fnorm));
}
