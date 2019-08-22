#include "bml.h"
#include "../macros.h"
#include "../typed.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>

static void TYPED_FUNC(
    ref_multiply) (
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
            C[ROWMAJOR(i, j, N, N)] *= beta;
            for (int k = 0; k < N; k++)
            {
                C[ROWMAJOR(i, j, N, N)] +=
                    alpha * A[ROWMAJOR(i, k, N, N)] * B[ROWMAJOR(k, j, N, N)];
            }
        }
    }
}

#if defined(SINGLE_REAL) || defined(SINGLE_COMPLEX)
#define ABS_TOL 2e-6
#else
#define ABS_TOL 1e-12
#endif

static int TYPED_FUNC(
    compare_matrix) (
    const int N,
    const bml_matrix_precision_t matrix_precision,
    const REAL_T * A,
    const REAL_T * B)
{
    for (int i = 0; i < N * N; i++)
    {
        if (ABS(A[i] - B[i]) > ABS_TOL)
        {
            bml_print_dense_matrix(N, matrix_precision, dense_row_major, A, 0,
                                   N, 0, N);
            bml_print_dense_matrix(N, matrix_precision, dense_row_major, B, 0,
                                   N, 0, N);
            LOG_INFO("element %d outside %1.2e\n", i, ABS_TOL);
#if defined(SINGLE_COMPLEX) || defined(DOUBLE_COMPLEX)
            LOG_INFO("A[%d] = %e+%ei\n", i, creal(A[i]), cimag(A[i]));
            LOG_INFO("B[%d] = %e+%ei\n", i, creal(B[i]), cimag(B[i]));
            LOG_INFO("abs(A[%d]-B[%d]) = %e\n", i, i, ABS(A[i] - B[i]));
#else
            LOG_INFO("A[%d] = %e\n", i, A[i]);
            LOG_INFO("B[%d] = %e\n", i, B[i]);
            LOG_INFO("abs(A[%d]-B[%d]) = %e\n", i, i, ABS(A[i] - B[i]));
#endif
            return 1;
        }
    }
    return 0;
}

int TYPED_FUNC(
    test_multiply) (
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

    A = bml_random_matrix(matrix_type, matrix_precision, N, M, sequential);
    B = bml_random_matrix(matrix_type, matrix_precision, N, M, sequential);
    C = bml_random_matrix(matrix_type, matrix_precision, N, M, sequential);

    bml_print_bml_matrix(A, 0, N, 0, N);
    bml_print_bml_matrix(B, 0, N, 0, N);
    bml_print_bml_matrix(C, 0, N, 0, N);

    A_dense = bml_export_to_dense(A, dense_row_major);
    B_dense = bml_export_to_dense(B, dense_row_major);
    C_dense = bml_export_to_dense(C, dense_row_major);
    D_dense = bml_export_to_dense(C, dense_row_major);

    bml_multiply(A, B, C, alpha, beta, threshold);
    E_dense = bml_export_to_dense(C, dense_row_major);

    TYPED_FUNC(ref_multiply) (N, A_dense, B_dense, D_dense, alpha, beta,
                              threshold);

    if (TYPED_FUNC(compare_matrix) (N, matrix_precision, D_dense, E_dense) !=
        0)
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
