#include "../../macros.h"
#include "../../typed.h"
#include "../bml_allocate.h"
#include "../bml_logger.h"
#include "../bml_parallel.h"
#include "../bml_submatrix.h"
#include "../bml_trace.h"
#include "../bml_types.h"
#include "bml_submatrix_ellsort.h"
#include "bml_trace_ellsort.h"
#include "bml_types_ellsort.h"

#include <complex.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/** Calculate the trace of a matrix.
 *
 *  \ingroup trace_group
 *
 *  \param A The matrix to calculate a trace for
 *  \return the trace of A
 */
double TYPED_FUNC(
    bml_trace_ellsort) (
    bml_matrix_ellsort_t * A)
{
    int N = A->N;
    int M = A->M;

    int *A_index = (int *) A->index;
    int *A_nnz = (int *) A->nnz;
    int *A_localRowMin = (int *) A->domain->localRowMin;
    int *A_localRowMax = (int *) A->domain->localRowMax;

    REAL_T trace = 0.0;
    REAL_T *A_value = (REAL_T *) A->value;

    int myRank = bml_getMyRank();

#pragma omp parallel for                        \
  shared(N, M, A_value, A_index, A_nnz)         \
  shared(A_localRowMin, A_localRowMax, myRank)  \
  reduction(+:trace)
    //for (int i = 0; i < N; i++)
    for (int i = A_localRowMin[myRank]; i < A_localRowMax[myRank]; i++)
    {
        for (int j = 0; j < A_nnz[i]; j++)
        {
            if (i == A_index[ROWMAJOR(i, j, N, M)])
            {
                trace += A_value[ROWMAJOR(i, j, N, M)];
                break;
            }
        }
    }

    return (double) REAL_PART(trace);
}

/** Calculate the trace of a matrix multiplication.
 * Both matrices must have the same size.
 *
 *  \ingroup trace_group
 *
 *  \param A The matrix A
 *  \param A The matrix B
 *  \return the trace of A*B
 */
double TYPED_FUNC(
    bml_trace_mult_ellsort) (
    bml_matrix_ellsort_t * A,
    bml_matrix_ellsort_t * B)
{
    int A_N = A->N;
    int A_M = A->M;

    int *A_index = (int *) A->index;
    int *A_nnz = (int *) A->nnz;
    int *A_localRowMin = (int *) A->domain->localRowMin;
    int *A_localRowMax = (int *) A->domain->localRowMax;

    REAL_T trace = 0.0;
    REAL_T *A_value = (REAL_T *) A->value;
    REAL_T *rvalue;

    int myRank = bml_getMyRank();

    if (A_N != B->N || A_M != B->M)
    {
        LOG_ERROR
            ("bml_trace_mult_ellsort: Matrices A and B have different sizes.");
    }

#pragma omp parallel for                        \
  private(rvalue)                               \
  shared(B, A_N, A_M, A_value, A_index, A_nnz)  \
  shared(A_localRowMin, A_localRowMax, myRank)  \
  reduction(+:trace)
    //for (int i = 0; i < A_N; i++)
    for (int i = A_localRowMin[myRank]; i < A_localRowMax[myRank]; i++)
    {
        rvalue =
            TYPED_FUNC(bml_getVector_ellsort) (B,
                                               &A_index[ROWMAJOR
                                                        (i, 0, A_N, A_M)], i,
                                               A_nnz[i]);
        for (int j = 0; j < A_nnz[i]; j++)
        {
            trace += A_value[ROWMAJOR(i, j, A_N, A_M)] * rvalue[j];
        }
        bml_free_memory(rvalue);
    }

    return (double) REAL_PART(trace);
}
