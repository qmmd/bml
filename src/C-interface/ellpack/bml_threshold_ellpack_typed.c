#include "../../macros.h"
#include "../../typed.h"
#include "../bml_allocate.h"
#include "../bml_parallel.h"
#include "../bml_threshold.h"
#include "../bml_types.h"
#include "bml_allocate_ellpack.h"
#include "bml_threshold_ellpack.h"
#include "bml_types_ellpack.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/** Threshold a matrix.
 *
 *  \ingroup threshold_group
 *
 *  \param A The matrix to be thresholded
 *  \param threshold Threshold value
 *  \return the thresholded A
 */
bml_matrix_ellpack_t
    * TYPED_FUNC(bml_threshold_new_ellpack) (bml_matrix_ellpack_t * A,
                                             double threshold)
{
    int N = A->N;
    int M = A->M;

    bml_matrix_ellpack_t *B =
        TYPED_FUNC(bml_zero_matrix_ellpack) (N, M, A->distribution_mode);

    REAL_T *A_value = (REAL_T *) A->value;
    int *A_index = A->index;
    int *A_nnz = A->nnz;
    int *A_localRowMin = A->domain->localRowMin;
    int *A_localRowMax = A->domain->localRowMax;

    REAL_T *B_value = (REAL_T *) B->value;
    int *B_index = B->index;
    int *B_nnz = B->nnz;

    int myRank = bml_getMyRank();

#pragma omp parallel for               \
    shared(N, M, A_value, A_index, A_nnz) \
    shared(A_localRowMin, A_localRowMax, myRank) \
    shared(B_value, B_index, B_nnz)
    //for (int i = 0; i < N; i++)
    for (int i = A_localRowMin[myRank]; i < A_localRowMax[myRank]; i++)
    {
        for (int j = 0; j < A_nnz[i]; j++)
        {
            if (is_above_threshold(A_value[ROWMAJOR(i, j, N, M)], threshold))
            {
                B_value[ROWMAJOR(i, B_nnz[i], N, M)] =
                    A_value[ROWMAJOR(i, j, N, M)];
                B_index[ROWMAJOR(i, B_nnz[i], N, M)] =
                    A_index[ROWMAJOR(i, j, N, M)];
                B_nnz[i]++;
            }
        }
    }

    return B;
}

/** Threshold a matrix in place.
 *
 *  \ingroup threshold_group
 *
 *  \param A The matrix to be thresholded
 *  \param threshold Threshold value
 *  \return the thresholded A
 */
void TYPED_FUNC(
    bml_threshold_ellpack) (
    bml_matrix_ellpack_t * A,
    double threshold)
{
    int N = A->N;
    int M = A->M;

    REAL_T *A_value = (REAL_T *) A->value;
    int *A_index = A->index;
    int *A_nnz = A->nnz;
    int *A_localRowMin = A->domain->localRowMin;
    int *A_localRowMax = A->domain->localRowMax;

    int myRank = bml_getMyRank();

    int rlen;
#pragma omp parallel for               \
    private(rlen) \
    shared(N,M,A_value,A_index,A_nnz) \
    shared(A_localRowMin, A_localRowMax, myRank)
    //for (int i = 0; i < N; i++)
    for (int i = A_localRowMin[myRank]; i < A_localRowMax[myRank]; i++)
    {
        rlen = 0;
        for (int j = 0; j < A_nnz[i]; j++)
        {
            if (is_above_threshold(A_value[ROWMAJOR(i, j, N, M)], threshold))
            {
                if (rlen < j)
                {
                    A_value[ROWMAJOR(i, rlen, N, M)] =
                        A_value[ROWMAJOR(i, j, N, M)];
                    A_index[ROWMAJOR(i, rlen, N, M)] =
                        A_index[ROWMAJOR(i, j, N, M)];
                }
                rlen++;
            }
        }
        A_nnz[i] = rlen;
    }
}
