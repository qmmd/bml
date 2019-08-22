#ifdef BML_USE_MAGMA
#include "magma_v2.h"
#include "../bml_allocate.h"
#endif

#include "../../macros.h"
#include "../../typed.h"
#include "../blas.h"
#include "../bml_norm.h"
#include "../bml_parallel.h"
#include "../bml_types.h"
#include "bml_norm_dense.h"
#include "bml_types_dense.h"

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
    bml_matrix_dense_t * A)
{
    int N = A->N;
    REAL_T sum = 0.0;

#ifdef BML_USE_MAGMA
#if defined(SINGLE_COMPLEX) || defined(DOUBLE_COMPLEX)
    MAGMA_T tsum = MAGMACOMPLEX(MAKE) (0., 0.);
    for (int i = 0; i < N; i++)
    {
        tsum =
            MAGMACOMPLEX(ADD) (tsum,
                               MAGMA(dotu) (N,
                                            (MAGMA_T *) A->matrix + i * A->ld,
                                            1,
                                            (MAGMA_T *) A->matrix + i * A->ld,
                                            1, A->queue));
    }
    sum = MAGMACOMPLEX(REAL) (tsum) + I * MAGMACOMPLEX(IMAG) (tsum);
#else
    for (int i = 0; i < N; i++)
    {
        sum += MAGMA(dot) (N, (MAGMA_T *) A->matrix + i * A->ld, 1,
                           (MAGMA_T *) A->matrix + i * A->ld, 1, A->queue);
    }
#endif

#else

    REAL_T *A_matrix = A->matrix;

    int *A_localRowMin = A->domain->localRowMin;
    int *A_localRowMax = A->domain->localRowMax;

    int myRank = bml_getMyRank();

#pragma omp parallel for                        \
  shared(N, A_matrix)                           \
  shared(A_localRowMin, A_localRowMax, myRank)  \
  reduction(+:sum)
    for (int i = A_localRowMin[myRank] * N; i < A_localRowMax[myRank] * N;
         i++)
    {
        sum += A_matrix[i] * A_matrix[i];
    }
#endif

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
    bml_matrix_dense_t * A,
    int core_size)
{
    int N = A->N;

    REAL_T sum = 0.0;

#ifdef BML_USE_MAGMA
#if defined(SINGLE_COMPLEX) || defined(DOUBLE_COMPLEX)
    MAGMA_T tsum = MAGMACOMPLEX(MAKE) (0., 0.);
    for (int i = 0; i < core_size; i++)
    {
        tsum =
            MAGMACOMPLEX(ADD) (tsum,
                               MAGMA(dotu) (N,
                                            (MAGMA_T *) A->matrix + i * A->ld,
                                            1,
                                            (MAGMA_T *) A->matrix + i * A->ld,
                                            1, A->queue));
    }
    sum = MAGMACOMPLEX(REAL) (tsum) + I * MAGMACOMPLEX(IMAG) (tsum);
#else
    for (int i = 0; i < core_size; i++)
    {
        sum += MAGMA(dot) (N, (MAGMA_T *) A->matrix + i * A->ld, 1,
                           (MAGMA_T *) A->matrix + i * A->ld, 1, A->queue);
    }
#endif

#else

    REAL_T *A_matrix = A->matrix;

#pragma omp parallel for                        \
  shared(N, A_matrix)                           \
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
#endif

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
    bml_matrix_dense_t * A,
    bml_matrix_dense_t * B,
    double alpha,
    double beta,
    double threshold)
{
    int N = A->N;

    REAL_T sum = 0.0;
#ifdef BML_USE_MAGMA            //do work on CPU for now...
    MAGMA_T *A_matrix = bml_allocate_memory(sizeof(MAGMA_T) * A->N * A->N);
    MAGMA(getmatrix) (A->N, A->N, A->matrix, A->ld, A_matrix, A->N, A->queue);
    MAGMA_T *B_matrix = bml_allocate_memory(sizeof(MAGMA_T) * B->N * B->N);
    MAGMA(getmatrix) (B->N, B->N, B->matrix, B->ld, B_matrix, B->N, B->queue);

#else
    REAL_T *A_matrix = A->matrix;
    REAL_T *B_matrix = B->matrix;
#endif

    int *A_localRowMin = A->domain->localRowMin;
    int *A_localRowMax = A->domain->localRowMax;

#ifdef BML_USE_MAGMA
    MAGMA_T alpha_ = MAGMACOMPLEX(MAKE) (alpha, 0.);
    MAGMA_T beta_ = MAGMACOMPLEX(MAKE) (beta, 0.);
#else
    REAL_T alpha_ = (REAL_T) alpha;
    REAL_T beta_ = (REAL_T) beta;
#endif

    int myRank = bml_getMyRank();

#pragma omp parallel for                        \
  shared(alpha_, beta_)                         \
  shared(N, A_matrix, B_matrix)                 \
  shared(A_localRowMin, A_localRowMax, myRank)  \
  reduction(+:sum)
    //for (int i = 0; i < N * N; i++)
    for (int i = A_localRowMin[myRank] * N; i < A_localRowMax[myRank] * N;
         i++)
    {
#ifdef BML_USE_MAGMA
        MAGMA_T ttemp =
            MAGMACOMPLEX(ADD) (MAGMACOMPLEX(MUL) (alpha_, A_matrix[i]),
                               MAGMACOMPLEX(MUL) (beta_, B_matrix[i]));
        REAL_T temp =
            MAGMACOMPLEX(REAL) (ttemp) + I * MAGMACOMPLEX(IMAG) (ttemp);
#else
        REAL_T temp = alpha_ * A_matrix[i] + beta_ * B_matrix[i];
#endif
        if (ABS(temp) > threshold)
            sum += temp * temp;
    }
#ifdef BML_USE_MAGMA
    free(A_matrix);
    free(B_matrix);
#endif
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
    bml_matrix_dense_t * A)
{
    double sum = 0.0;
//#ifdef BML_USE_MAGMA
//    int lwork=A->N;
//    MAGMA_T *dwork;
//    MAGMA(malloc)(&dwork,lwork);
//    sum = MAGMABLAS(lange)(MagmaFrobeniusNorm, A->N, A->N, A->matrix, A->ld,
//                           dwork, lwork, A->queue);
//#else
    sum = TYPED_FUNC(bml_sum_squares_dense) (A);

#ifdef DO_MPI
    if (bml_getNRanks() > 1 && A->distribution_mode == distributed)
    {
        bml_sumRealReduce(&sum);
    }
#endif
//#endif

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
    bml_matrix_dense_t * A,
    bml_matrix_dense_t * B)
{
    int N = A->N;

    REAL_T *A_matrix = (REAL_T *) A->matrix;
    REAL_T *B_matrix = (REAL_T *) B->matrix;

    int *A_localRowMin = A->domain->localRowMin;
    int *A_localRowMax = A->domain->localRowMax;

    double fnorm = 0.0;
    REAL_T temp;

    int myRank = bml_getMyRank();

#pragma omp parallel for                        \
  shared(temp)                                  \
  shared(N, A_matrix, B_matrix)                 \
  shared(A_localRowMin, A_localRowMax, myRank)  \
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
