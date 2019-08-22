#include "../../macros.h"
#include "../../typed.h"
#include "../bml_allocate.h"
#include "../bml_copy.h"
#include "../bml_types.h"
#include "bml_allocate_ellpack.h"
#include "bml_copy_ellpack.h"
#include "bml_types_ellpack.h"

#include <complex.h>
#include <stdlib.h>
#include <string.h>

/** Copy an ellpack matrix - result is a new matrix.
 *
 *  \ingroup copy_group
 *
 *  \param A The matrix to be copied
 *  \return A copy of matrix A.
 */
bml_matrix_ellpack_t *TYPED_FUNC(
    bml_copy_ellpack_new) (
    bml_matrix_ellpack_t * A)
{
    bml_matrix_dimension_t matrix_dimension = { A->N, A->N, A->M };
    bml_matrix_ellpack_t *B =
        TYPED_FUNC(bml_noinit_matrix_ellpack) (matrix_dimension,
                                               A->distribution_mode);

    int N = A->N;
    int M = A->M;

    int *A_index = A->index;
    int *A_nnz = A->nnz;
    REAL_T *A_value = A->value;

    int *B_index = B->index;
    int *B_nnz = B->nnz;
    REAL_T *B_value = B->value;

#ifdef NOGPU
    // copy to host and copy back
#pragma omp target update from(A_nnz[:N], A_index[:N*M], A_value[:N*M])
    //    memcpy(B->index, A->index, sizeof(int) * A->N * A->M);
    memcpy(B->nnz, A->nnz, sizeof(int) * A->N);
    //    memcpy(B->value, A->value, sizeof(REAL_T) * A->N * A->M);
#pragma omp parallel for
    for (int i = 0; i < N; i++)
    {
        memcpy(&B_index[ROWMAJOR(i, 0, N, M)], &A_index[ROWMAJOR(i, 0, N, M)],
               M * sizeof(int));
        memcpy(&B_value[ROWMAJOR(i, 0, N, M)], &A_value[ROWMAJOR(i, 0, N, M)],
               M * sizeof(REAL_T));
        //      A_nnz[perm[i]] = B_nnz[i];
    }
// push the data to the GPU
#pragma omp target update to(B_nnz[:N], B_index[:N*M], B_value[:N*M])
#endif // NOGPU

    // all work on device

#pragma omp target teams distribute parallel for 
    for (int i = 0; i < N; i++)
    {
        B_nnz[i] = A_nnz[i];
    }

#pragma omp target teams distribute parallel for collapse(2) schedule (static, 1)
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
            B_index[ROWMAJOR(i,j,N,M)] = A_index[ROWMAJOR(i, j, N, M)];
            B_value[ROWMAJOR(i,j,N,M)] = A_value[ROWMAJOR(i, j, N, M)];
        }
    }

    bml_copy_domain(A->domain, B->domain);
    bml_copy_domain(A->domain2, B->domain2);
    return B;
}

/** Copy an ellpack matrix.
 *
 *  \ingroup copy_group
 *
 *  \param A The matrix to be copied
 *  \param B Copy of matrix A
 */
void TYPED_FUNC(
    bml_copy_ellpack) (
    bml_matrix_ellpack_t * A,
    bml_matrix_ellpack_t * B)
{
    int N = A->N;
    int M = A->M;

    int *A_index = A->index;
    int *A_nnz = A->nnz;
    REAL_T *A_value = A->value;

    int *B_index = B->index;
    int *B_nnz = B->nnz;
    REAL_T *B_value = B->value;

#ifdef NOGPU
#pragma omp target update from(A_nnz[:N], A_index[:N*M], A_value[:N*M])
    // memcpy(B->index, A->index, sizeof(int) * A->N * A->M);
    memcpy(B->nnz, A->nnz, sizeof(int) * A->N);
    //    memcpy(B->value, A->value, sizeof(REAL_T) * A->N * A->M);
#pragma omp parallel for
    for (int i = 0; i < N; i++)
    {
        memcpy(&B_index[ROWMAJOR(i, 0, N, M)], &A_index[ROWMAJOR(i, 0, N, M)],
               M * sizeof(int));
        memcpy(&B_value[ROWMAJOR(i, 0, N, M)], &A_value[ROWMAJOR(i, 0, N, M)],
               M * sizeof(REAL_T));
        //      A_nnz[perm[i]] = B_nnz[i];
    }
// push the data to the GPU
#pragma omp target update to(B_nnz[:N], B_index[:N*M], B_value[:N*M])

#endif // NOGPU

// All data and copy stays on deveice
#pragma omp target teams distribute parallel for 
    for (int i = 0; i < N; i++)
    {
        B_nnz[i] = A_nnz[i];
    }

#pragma omp target teams distribute parallel for collapse(2) schedule (static, 1)
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
            B_index[ROWMAJOR(i,j,N,M)] = A_index[ROWMAJOR(i, j, N, M)];
            B_value[ROWMAJOR(i,j,N,M)] = A_value[ROWMAJOR(i, j, N, M)];
        }
    }


    if (A->distribution_mode == B->distribution_mode)
    {
        bml_copy_domain(A->domain, B->domain);
        bml_copy_domain(A->domain2, B->domain2);
    }
}

/** Reorder an ellpack matrix.
 *
 *  \ingroup copy_group
 *
 *  \param A The matrix to be reordered
 *  \param B The permutation vector
 */
void TYPED_FUNC(
    bml_reorder_ellpack) (
    bml_matrix_ellpack_t * A,
    int *perm)
{
    int N = A->N;
    int M = A->M;

    int *A_index = A->index;
    int *A_nnz = A->nnz;
    REAL_T *A_value = A->value;

    bml_matrix_ellpack_t *B = bml_copy_new(A);
    int *B_index = B->index;
    int *B_nnz = B->nnz;
    REAL_T *B_value = B->value;

    // Reorder rows - need to copy
#pragma omp parallel for
    for (int i = 0; i < N; i++)
    {
        memcpy(&A_index[ROWMAJOR(perm[i], 0, N, M)],
               &B_index[ROWMAJOR(i, 0, N, M)], M * sizeof(int));
        memcpy(&A_value[ROWMAJOR(perm[i], 0, N, M)],
               &B_value[ROWMAJOR(i, 0, N, M)], M * sizeof(REAL_T));
        A_nnz[perm[i]] = B_nnz[i];
    }

    bml_deallocate_ellpack(B);

    // Reorder elements in each row - just change index
#pragma omp parallel for
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < A_nnz[i]; j++)
        {
            A_index[ROWMAJOR(i, j, N, M)] =
                perm[A_index[ROWMAJOR(i, j, N, M)]];
        }
    }
}
