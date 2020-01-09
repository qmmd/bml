#include "bml.h"
#include "bml_test.h"

#include "../typed.h"
#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

int TYPED_FUNC(
    test_io_matrix) (
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
    bml_write_bml_matrix(A, "ctest_matrix.mtx");
    B = bml_zero_matrix(matrix_type, matrix_precision, N, M, sequential);
    bml_read_bml_matrix(B, "ctest_matrix.mtx");

    A_dense = bml_export_to_dense(A, dense_row_major);
    B_dense = bml_export_to_dense(B, dense_row_major);
    bml_print_dense_matrix(N, matrix_precision, dense_row_major, A_dense, 0,
                           N, 0, N);
    bml_print_dense_matrix(N, matrix_precision, dense_row_major, B_dense, 0,
                           N, 0, N);
    for (int i = 0; i < N * N; i++)
    {
        if (fabs(A_dense[i] - B_dense[i]) > 1e-12)
        {
            LOG_ERROR("matrices are not identical; A[%d] = %e\n", i,
                      A_dense[i]);
            return -1;
        }
    }
    bml_free_memory(A_dense);
    bml_free_memory(B_dense);
    bml_deallocate(&A);
    bml_deallocate(&B);

    LOG_INFO("io matrix test passed\n");

    return 0;
}
