#include "../../macros.h"
#include "../../typed.h"
#include "bml_allocate.h"
#include "bml_allocate_dense.h"
#include "bml_transpose.h"
#include "bml_parallel.h"
#include "bml_transpose_dense.h"
#include "bml_types.h"
#include "bml_types_dense.h"

#include <complex.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/** Transpose a matrix.
 *
 *  \ingroup transpose_group
 *
 *  \param A The matrix to be transposed
 *  \return The transposed A
 */
bml_matrix_dense_t *TYPED_FUNC(
    bml_transpose_new_dense) (
    const bml_matrix_dense_t * A)
{
    int N = A->N;

    bml_matrix_dense_t *B =
        TYPED_FUNC(bml_zero_matrix_dense) (N, A->distribution_mode);
    REAL_T *A_matrix = A->matrix;
    REAL_T *B_matrix = B->matrix;

    int *A_localRowMin = A->domain->localRowMin;
    int *A_localRowMax = A->domain->localRowMax;

    int myRank = bml_getMyRank();

#pragma omp parallel for default(none) \
    shared(N, A_matrix, B_matrix) \
    shared(A_localRowMin, A_localRowMax, myRank)
    //for (int i = 0; i < N; i++)
    for (int i = A_localRowMin[myRank]; i < A_localRowMax[myRank]; i++)
    {
        for (int j = 0; j < N; j++)
        {
            B_matrix[LINEARINDEX(i, j, N, N)] =
                A_matrix[LINEARINDEX(j, i, N, N)];
        }
    }
    return B;
}

/** Transpose a matrix in place.
 *
 *  \ingroup transpose_group
 *
 *  \param A The matrix to be transposed
 *  \return The transposed A
 */
void TYPED_FUNC(
    bml_transpose_dense) (
    bml_matrix_dense_t * A)
{
    int N = A->N;

    REAL_T *A_matrix = A->matrix;
    REAL_T tmp;

#pragma omp parallel for default(none) \
    private(tmp) \
    shared(N, A_matrix)
    for (int i = 0; i < N - 1; i++)
    {
        for (int j = i + 1; j < N; j++)
        {
            if (i != j)
            {
                tmp = A_matrix[LINEARINDEX(i, j, N, N)];
                A_matrix[LINEARINDEX(i, j, N, N)] =
                    A_matrix[LINEARINDEX(j, i, N, N)];
                A_matrix[LINEARINDEX(j, i, N, N)] = tmp;
            }
        }
    }
}
