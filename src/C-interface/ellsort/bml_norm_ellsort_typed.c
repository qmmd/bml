#include "../../macros.h"
#include "../../typed.h"
#include "../bml_norm.h"
#include "../bml_parallel.h"
#include "../bml_types.h"
#include "bml_norm_ellsort.h"
#include "bml_types_ellsort.h"

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
    bml_sum_squares_ellsort) (
    bml_matrix_ellsort_t * A)
{
    int N = A->N;
    int M = A->M;

    int *A_nnz = (int *) A->nnz;
    int *A_localRowMin = A->domain->localRowMin;
    int *A_localRowMax = A->domain->localRowMax;

    REAL_T sum = 0.0;
    REAL_T *A_value = (REAL_T *) A->value;

    int myRank = bml_getMyRank();

#pragma omp parallel for                        \
  shared(N, M, A_value, A_nnz)                  \
  shared(A_localRowMin, A_localRowMax, myRank)  \
  reduction(+:sum)

    //for (int i = 0; i < N; i++)
    for (int i = A_localRowMin[myRank]; i < A_localRowMax[myRank]; i++)
    {
        for (int j = 0; j < A_nnz[i]; j++)
        {
            REAL_T xval = A_value[ROWMAJOR(i, j, N, M)];
            sum += xval * xval;
        }
    }

    return (double) REAL_PART(sum);
}

/** Calculate the sum of squares of all the core elements of a submatrix.
 *
 *  \ingroup norm_group
 *
 *  \param A The matrix
 *  \param core_pos Core rows of submatrix
 *  \param core_size Number of core rows
 *  \return The sum of squares of A
 */
double TYPED_FUNC(
    bml_sum_squares_submatrix_ellsort) (
    bml_matrix_ellsort_t * A,
    int core_size)
{
    int N = A->N;
    int M = A->M;

    int *A_index = (int *) A->index;
    int *A_nnz = (int *) A->nnz;

    REAL_T sum = 0.0;
    REAL_T *A_value = (REAL_T *) A->value;

#pragma omp parallel for                        \
  shared(N, M, A_index, A_nnz, A_value)         \
  reduction(+:sum)
    for (int i = 0; i < core_size; i++)
    {
        for (int j = 0; j < A_nnz[i]; j++)
        {
            if (A_index[ROWMAJOR(i, j, N, M)] < core_size)
            {
                REAL_T value = A_value[ROWMAJOR(i, j, N, M)];
                sum += value * value;
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
    bml_sum_squares2_ellsort) (
    bml_matrix_ellsort_t * A,
    bml_matrix_ellsort_t * B,
    double alpha,
    double beta,
    double threshold)
{
    int A_N = A->N;
    int A_M = A->M;
    int B_N = B->N;
    int B_M = B->M;

    int *A_index = (int *) A->index;
    int *A_nnz = (int *) A->nnz;
    int *B_index = (int *) B->index;
    int *B_nnz = (int *) B->nnz;

    int *A_localRowMin = A->domain->localRowMin;
    int *A_localRowMax = A->domain->localRowMax;

    REAL_T sum = 0.0;
    REAL_T *A_value = (REAL_T *) A->value;
    REAL_T *B_value = (REAL_T *) B->value;

    REAL_T alpha_ = (REAL_T) alpha;
    REAL_T beta_ = (REAL_T) beta;

    int myRank = bml_getMyRank();

#if !(defined(__IBMC__) || defined(__ibmxl__))
    REAL_T y[A_N];
    int ix[A_N], jjb[A_N];

    memset(y, 0.0, A_N * sizeof(REAL_T));
    memset(ix, 0, A_N * sizeof(int));
    memset(jjb, 0, A_N * sizeof(int));
#endif

#if defined(__IBMC__) || defined(__ibmxl__)
#pragma omp parallel for                        \
  shared(alpha_, beta_)                         \
  shared(A_N, A_M, A_index, A_nnz, A_value)     \
  shared(A_localRowMin, A_localRowMax, myRank)  \
  shared(B_N, B_M, B_index, B_nnz, B_value)     \
  reduction(+:sum)
#else
#pragma omp parallel for                        \
  shared(alpha_, beta_)                         \
  shared(A_N, A_M, A_index, A_nnz, A_value)     \
  shared(A_localRowMin, A_localRowMax, myRank)  \
  shared(B_N, B_M, B_index, B_nnz, B_value)     \
  firstprivate(ix, jjb, y) \
  reduction(+:sum)
#endif

    //for (int i = 0; i < A_N; i++)
    for (int i = A_localRowMin[myRank]; i < A_localRowMax[myRank]; i++)
    {

#if defined(__IBMC__) || defined(__ibmxl__)
        REAL_T y[A_N];
        int ix[A_N], jjb[A_N];

        memset(ix, 0, A_N * sizeof(int));
#endif

        int l = 0;
        for (int jp = 0; jp < A_nnz[i]; jp++)
        {
            int k = A_index[ROWMAJOR(i, jp, A_N, A_M)];
            if (ix[k] == 0)
            {
                y[k] = 0.0;
                ix[k] = i + 1;
                jjb[l] = k;
                l++;
            }
            y[k] += alpha_ * A_value[ROWMAJOR(i, jp, A_N, A_M)];
        }

        for (int jp = 0; jp < B_nnz[i]; jp++)
        {
            int k = B_index[ROWMAJOR(i, jp, B_N, B_M)];
            if (ix[k] == 0)
            {
                y[k] = 0.0;
                ix[k] = i + 1;
                jjb[l] = k;
                l++;
            }
            y[k] += beta_ * B_value[ROWMAJOR(i, jp, B_N, B_M)];
        }

        for (int jp = 0; jp < l; jp++)
        {
            if (ABS(y[jjb[jp]]) > threshold)
                sum += y[jjb[jp]] * y[jjb[jp]];

            ix[jjb[jp]] = 0;
            y[jjb[jp]] = 0.0;
            jjb[jp] = 0;
        }
    }

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
    bml_fnorm_ellsort) (
    bml_matrix_ellsort_t * A)
{
    double fnorm = TYPED_FUNC(bml_sum_squares_ellsort) (A);
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
    bml_fnorm2_ellsort) (
    bml_matrix_ellsort_t * A,
    bml_matrix_ellsort_t * B)
{
    int N = A->N;
    int M = A->M;
    double fnorm = 0.0;
    REAL_T rvalue;

    int *A_nnz = (int *) A->nnz;
    int *A_index = (int *) A->index;
    int *A_localRowMin = A->domain->localRowMin;
    int *A_localRowMax = A->domain->localRowMax;
    REAL_T *A_value = (REAL_T *) A->value;
    int *B_nnz = (int *) B->nnz;
    int *B_index = (int *) B->index;
    REAL_T *B_value = (REAL_T *) B->value;

    REAL_T temp;

    int myRank = bml_getMyRank();

#pragma omp parallel for                        \
  private(rvalue, temp)                         \
  shared(N, M, A_nnz, A_index, A_value)         \
  shared(A_localRowMin, A_localRowMax, myRank)  \
  shared(B_nnz, B_index, B_value)               \
  reduction(+:fnorm)
    //for (int i = 0; i < N; i++)
    for (int i = A_localRowMin[myRank]; i < A_localRowMax[myRank]; i++)
    {
        for (int j = 0; j < A_nnz[i]; j++)
        {
            for (int k = 0; k < B_nnz[i]; k++)
            {
                if (A_index[ROWMAJOR(i, j, N, M)] ==
                    B_index[ROWMAJOR(i, k, N, M)])
                {
                    rvalue = B_value[ROWMAJOR(i, k, N, M)];
                    break;
                }
                rvalue = 0.0;
            }

            temp = A_value[ROWMAJOR(i, j, N, M)] - rvalue;
            fnorm += temp * temp;
        }

        for (int j = 0; j < B_nnz[i]; j++)
        {
            for (int k = 0; k < A_nnz[i]; k++)
            {
                if (A_index[ROWMAJOR(i, k, N, M)] ==
                    B_index[ROWMAJOR(i, j, N, M)])
                {
                    rvalue = A_value[ROWMAJOR(i, k, N, M)];
                    break;
                }
                rvalue = 0.0;
            }

            if (rvalue == 0.0)
            {
                temp = B_value[ROWMAJOR(i, j, N, M)];
                fnorm += temp * temp;
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

    return (double) REAL_PART(fnorm);
}
