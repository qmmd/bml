#include "bml.h"
#include "bml_test.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>

void
multiply(
    const int N,
    const REAL_T * A,
    const REAL_T * B,
    REAL_T * C,
    const double alpha,
    const double beta,
    const double threshold)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            C[i * N + j] *= beta;
            for (int k = 0; k < N; k++)
            {
                C[i * N + j] += alpha * A[i * N + k] * B[k * N + j];
            }
        }
    }
}

#if defined(SINGLE_REAL) || defined(SINGLE_COMPLEX)
#define ABS_TOL 1e-6
#else
#define ABS_TOL 1e-12
#endif

int
compare_matrix(
    const int N,
    const bml_matrix_precision_t matrix_precision,
    const REAL_T * A,
    const REAL_T * B)
{
    for (int i = 0; i < N * N; i++)
    {
        if (fabs(A[i] - B[i]) > ABS_TOL)
        {
            bml_print_dense_matrix(N, matrix_precision, dense_row_major, A, 0,
                                   N, 0, N);
            bml_print_dense_matrix(N, matrix_precision, dense_row_major, B, 0,
                                   N, 0, N);
            return 1;
        }
    }
    return 0;
}

int
test_function(
    const int N,
    const bml_matrix_type_t matrix_type,
    const bml_matrix_precision_t matrix_precision,
    const int M)
{
    bml_matrix_t *A = NULL;
    bml_matrix_t *B = NULL;
    bml_matrix_t *C = NULL;

    REAL_T *A_dense = NULL;
    REAL_T *B_dense = NULL;
    REAL_T *C_dense = NULL;
    REAL_T *D_dense = NULL;
    REAL_T *E_dense = NULL;

    const double alpha = 1.2;
    const double beta = 0.8;
    const double threshold = 0.0;

    A = bml_banded_matrix(matrix_type, matrix_precision, N, 3);
    B = bml_banded_matrix(matrix_type, matrix_precision, N, 3);
    C = bml_banded_matrix(matrix_type, matrix_precision, N, M);

    bml_print_bml_matrix(A, 0, N, 0, N);
    bml_print_bml_matrix(B, 0, N, 0, N);
    bml_print_bml_matrix(C, 0, N, 0, N);

    A_dense = bml_export_to_dense(A, dense_row_major);
    B_dense = bml_export_to_dense(B, dense_row_major);
    C_dense = bml_export_to_dense(C, dense_row_major);
    D_dense = bml_export_to_dense(C, dense_row_major);

    bml_multiply(A, B, C, alpha, beta, threshold);
    E_dense = bml_export_to_dense(C, dense_row_major);

    multiply(N, A_dense, B_dense, D_dense, alpha, beta, threshold);

    if (compare_matrix(N, matrix_precision, D_dense, E_dense) != 0)
    {
        LOG_ERROR("matrix product incorrect\n");
        return -1;
    }
    LOG_INFO("multiply matrix test passed\n");

    bml_deallocate(&A);
    bml_deallocate(&B);
    bml_deallocate(&C);

    bml_free_memory(A_dense);
    bml_free_memory(B_dense);
    bml_free_memory(C_dense);
    bml_free_memory(D_dense);
    bml_free_memory(E_dense);

    return 0;
}
