#include "../../typed.h"
#include "../bml_allocate.h"
#include "../bml_threshold.h"
#include "../bml_parallel.h"
#include "../bml_types.h"
#include "bml_allocate_dense.h"
#include "bml_threshold_dense.h"
#include "bml_types_dense.h"
#include "../bml_logger.h"

#include <complex.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/** Threshold a matrix.
 *
 *  \ingroup threshold_group
 *
 *  \param A The matrix to be thresholded
 *  \param threshold Threshold value
 *  \return The thresholded A
 */
bml_matrix_dense_t *TYPED_FUNC(
    bml_threshold_new_dense) (
    bml_matrix_dense_t * A,
    double threshold)
{
#ifdef BML_USE_MAGMA
    LOG_ERROR
        ("bml_threshold_new_dense() not implemented for MAGMA matrices\n");
#endif

    int N = A->N;
    bml_matrix_dimension_t matrix_dimension = { A->N, A->N, A->N };
    bml_matrix_dense_t *B =
        TYPED_FUNC(bml_zero_matrix_dense) (matrix_dimension,
                                           A->distribution_mode);
    REAL_T *A_matrix = A->matrix;
    REAL_T *B_matrix = B->matrix;

    int *A_localRowMin = A->domain->localRowMin;
    int *A_localRowMax = A->domain->localRowMax;

    int myRank = bml_getMyRank();

#pragma omp parallel for                        \
  shared(N, A_matrix, B_matrix)                 \
  shared(A_localRowMin, A_localRowMax, myRank)
    //for (int i = 0; i < N * N; i++)
    for (int i = A_localRowMin[myRank] * N; i < A_localRowMax[myRank] * N;
         i++)
    {
        if (is_above_threshold(A_matrix[i], (REAL_T) threshold))
        {
            B_matrix[i] = A_matrix[i];
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
 *  \return The thresholded A
 */
void TYPED_FUNC(
    bml_threshold_dense) (
    bml_matrix_dense_t * A,
    double threshold)
{
#ifdef BML_USE_MAGMA
    LOG_ERROR("bml_threshold_dense() not implemented for MAGMA matrices\n");
#endif

    int N = A->N;
    REAL_T *A_matrix = A->matrix;

    int *A_localRowMin = A->domain->localRowMin;
    int *A_localRowMax = A->domain->localRowMax;

    int myRank = bml_getMyRank();

#pragma omp parallel for                        \
  shared(N, A_matrix)                           \
  shared(A_localRowMin, A_localRowMax, myRank)
    //for (int i = 0; i < N * N; i++)
    for (int i = A_localRowMin[myRank] * N; i < A_localRowMax[myRank] * N;
         i++)
    {
        if (!is_above_threshold(A_matrix[i], (REAL_T) threshold))
        {
            A_matrix[i] = (REAL_T) 0.0;
        }
    }
}
