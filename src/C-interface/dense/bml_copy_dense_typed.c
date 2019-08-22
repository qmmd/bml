#include "../../macros.h"
#include "../../typed.h"
#include "../bml_allocate.h"
#include "../bml_copy.h"
#include "../bml_types.h"
#include "bml_allocate_dense.h"
#include "bml_copy_dense.h"
#include "bml_types_dense.h"

#ifdef BML_USE_MAGMA
#include "magma_v2.h"
#endif

#include <complex.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/** Copy a dense matrix - result in new matrix.
 *
 *  \ingroup copy_group
 *
 *  \param A The matrix to be copied
 *  \return A copy of matrix A.
 */
bml_matrix_dense_t *TYPED_FUNC(
    bml_copy_dense_new) (
    bml_matrix_dense_t * A)
{
    bml_matrix_dimension_t matrix_dimension = { A->N, A->N, A->N };
    bml_matrix_dense_t *B =
        TYPED_FUNC(bml_zero_matrix_dense) (matrix_dimension,
                                           A->distribution_mode);
#ifdef BML_USE_MAGMA
    MAGMA(copymatrix) (A->N, A->N, A->matrix, A->ld,
                       B->matrix, B->ld, A->queue);
#else
    memcpy(B->matrix, A->matrix, sizeof(REAL_T) * A->N * A->N);
#endif
    bml_copy_domain(A->domain, B->domain);
    bml_copy_domain(A->domain2, B->domain2);
    return B;
}

/** Copy a dense matrix.
 *
 *  \ingroup copy_group
 *
 *  \param A The matrix to be copied
 *  \param B Copy of matrix A
 */
void TYPED_FUNC(
    bml_copy_dense) (
    bml_matrix_dense_t * A,
    bml_matrix_dense_t * B)
{
#ifdef BML_USE_MAGMA
    MAGMA(copymatrix) (A->N, A->N, A->matrix, A->ld,
                       B->matrix, B->ld, A->queue);
#else
    memcpy(B->matrix, A->matrix, sizeof(REAL_T) * A->N * A->N);
#endif
    if (A->distribution_mode == B->distribution_mode)
    {
        bml_copy_domain(A->domain, B->domain);
        bml_copy_domain(A->domain2, B->domain2);
    }
}

/** Reorder a dense matrix using a permutation vector.
 *
 *  \ingroup copy_group
 *
 *  \param A The matrix to be reordered
 *  \param perm The permutation vector
 */
void TYPED_FUNC(
    bml_reorder_dense) (
    bml_matrix_dense_t * A,
    int *perm)
{
    int N = A->N;

    bml_matrix_dense_t *B = TYPED_FUNC(bml_copy_dense_new) (A);
    REAL_T *A_matrix = A->matrix;
    REAL_T *B_matrix = B->matrix;

    // Reorder rows - need to copy
#pragma omp parallel for
    for (int i = 0; i < N; i++)
    {
        memcpy(&A_matrix[ROWMAJOR(perm[i], 0, N, N)],
               &B_matrix[ROWMAJOR(i, 0, N, N)], N * sizeof(REAL_T));
    }

    // Reorder elements in each row - just change index
#pragma omp parallel for
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            A_matrix[ROWMAJOR(i, j, N, N)] =
                B_matrix[ROWMAJOR(i, perm[j], N, N)];
        }
    }

    bml_deallocate_dense(B);
}
