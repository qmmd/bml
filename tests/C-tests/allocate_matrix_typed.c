#include "bml.h"
#include "../typed.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>

int TYPED_FUNC(
    test_allocate) (
    const int N,
    const bml_matrix_type_t matrix_type,
    const bml_matrix_precision_t matrix_precision,
    const int M)
{
    bml_matrix_t *A = NULL;
    bml_matrix_t *B = NULL;
    REAL_T *A_dense = NULL;
    REAL_T *B_dense = NULL;

    A = bml_random_matrix(matrix_type, matrix_precision, N, M, sequential);
    A_dense = bml_export_to_dense(A, dense_row_major);
    B = bml_import_from_dense(matrix_type, matrix_precision, dense_row_major,
                              N, M, A_dense, 0, sequential);
    B_dense = bml_export_to_dense(B, dense_row_major);
    bml_print_dense_matrix(N, matrix_precision, dense_row_major, A_dense, 0,
                           N, 0, N);
    bml_print_dense_matrix(N, matrix_precision, dense_row_major, B_dense, 0,
                           N, 0, N);
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if (i == j)
            {
                if (ABS(A_dense[i * N + j] - B_dense[i * N + j]) > 1e-12)
                {
                    LOG_ERROR
                        ("incorrect value on diagonal; A[%d,%d] = %e B[%d,%d] = %e\n",
                         i, i, A_dense[i * N + j], i, i, B_dense[i * N + j]);
                    return -1;
                }
            }
            else
            {
                if (ABS(A_dense[i * N + j] - B_dense[i * N + j]) > 1e-12)
                {
                    LOG_ERROR
                        ("incorrect value off-diagonal; A[%d,%d] = %e B[%d,%d] = %e\n",
                         i, j, A_dense[i * N + j], i, i, B_dense[i * N + j]);
                    return -1;
                }
            }
        }
    }
    LOG_INFO("random matrix test passed\n");
    bml_free_memory(A_dense);
    bml_free_memory(B_dense);
    bml_deallocate(&A);
    bml_deallocate(&B);
    A = bml_identity_matrix(matrix_type, matrix_precision, N, M, sequential);
    A_dense = bml_export_to_dense(A, dense_row_major);
    bml_print_dense_matrix(N, matrix_precision, dense_row_major, A_dense, 0,
                           N, 0, N);
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if (i == j)
            {
                if (ABS(A_dense[i * N + j] - 1) > 1e-12)
                {
                    LOG_ERROR
                        ("incorrect value on diagonal; A[%d,%d] = %e\n",
                         i, i, A_dense[i * N + j]);
                    return -1;
                }
            }
            else
            {
                if (ABS(A_dense[i * N + j]) > 1e-12)
                {
                    LOG_ERROR
                        ("incorrect value off-diagonal; A[%d,%d] = %e\n",
                         i, j, A_dense[i * N + j]);
                    return -1;
                }
            }
        }
    }
    bml_free_memory(A_dense);
    LOG_INFO("identity matrix test passed\n");

    bml_clear(A);
    bml_deallocate(&A);

    return 0;
}
