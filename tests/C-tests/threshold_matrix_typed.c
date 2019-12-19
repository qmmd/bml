#include "bml.h"
#include "../typed.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

int TYPED_FUNC(
    test_threshold) (
    const int N,
    const bml_matrix_type_t matrix_type,
    const bml_matrix_precision_t matrix_precision,
    const int M)
{
    bml_matrix_t *A;
    bml_matrix_t *B;

    REAL_T *A_dense = NULL;
    REAL_T *B_dense = NULL;

    double threshold = 0.5;

#ifdef BML_USE_MAGMA
    if (matrix_type == dense)
    {
        LOG_INFO("submatrix matrix test not available\n");
        return 0;
    }
#endif

    A = bml_random_matrix(matrix_type, matrix_precision, N, M, sequential);
    B = bml_threshold_new(A, threshold);

    A_dense = bml_export_to_dense(A, dense_row_major);

    B_dense = bml_export_to_dense(B, dense_row_major);

    printf("A = \n");
    bml_print_dense_matrix(N, matrix_precision, dense_row_major, A_dense, 0,
                           N, 0, N);
    printf("B = threshold_new(A, %f)\n", threshold);
    bml_print_dense_matrix(N, matrix_precision, dense_row_major, B_dense, 0,
                           N, 0, N);
    for (int i = 0; i < N * N; i++)
    {
        if (ABS(B_dense[i]) > 0 && ABS(B_dense[i]) < threshold)
        {
            LOG_ERROR("matrices not thresholded B[%d] = %e\n", i, B_dense[i]);
            return -1;
        }
    }
    bml_free_memory(A_dense);
    bml_free_memory(B_dense);
    bml_deallocate(&A);
    bml_deallocate(&B);

    LOG_INFO("threshold matrix test passed\n");

    return 0;
}
