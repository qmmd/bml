#include "bml.h"
#include "../internal-blas/bml_gemm.h"
#include "../macros.h"
#include "../typed.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

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
    for (int j = 0; j < N; j++)
    {
        for (int i = 0; i < N; i++)
        {
            C[COLMAJOR(i, j, N, N)] *= beta;
            for (int k = 0; k < N; k++)
            {
                C[COLMAJOR(i, j, N, N)] +=
                    alpha * A[COLMAJOR(i, k, N, N)] * B[COLMAJOR(k, j, N, N)];
            }
        }
    }
}

/** B <- transpose(A)
 */
static void TYPED_FUNC(
    transpose) (
    const int N,
    const REAL_T * A,
    REAL_T * B)
{
    for (int j = 0; j < N; j++)
    {
        for (int i = 0; i < N; i++)
        {
            B[COLMAJOR(i, j, N, N)] = A[COLMAJOR(j, i, N, N)];
        }
    }
}

static void TYPED_FUNC(
    ref_random_matrix) (
    const int N,
    REAL_T * A)
{
    for (int j = 0; j < N; j++)
    {
        for (int i = 0; i < N; i++)
        {
            A[COLMAJOR(i, j, N, N)] =
                (REAL_T) (rand() / (double) RAND_MAX - 0.5);
        }
    }
}

/** Copy matrix A -> B.
 */
static void TYPED_FUNC(
    copy_matrix) (
    const int N,
    const REAL_T * A,
    REAL_T * B)
{
    memcpy(B, A, sizeof(REAL_T) * N * N);
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
            LOG_INFO("compare A:\n");
            bml_print_dense_matrix(N, matrix_precision, dense_row_major, A, 0,
                                   N, 0, N);
            LOG_INFO("compare B:\n");
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
    test_bml_gemm) (
    const int N,
    const bml_matrix_type_t matrix_type,
    const bml_matrix_precision_t matrix_precision,
    const int M)
{
    REAL_T alpha = 1.2;
    REAL_T beta = 0.8;

    char *transa[] = { "N", "T" };
    char *transb[] = { "N", "T" };

    REAL_T *A = calloc(sizeof(REAL_T), N * N);
    REAL_T *A_input = calloc(sizeof(REAL_T), N * N);
    REAL_T *B = calloc(sizeof(REAL_T), N * N);
    REAL_T *B_input = calloc(sizeof(REAL_T), N * N);
    REAL_T *C = calloc(sizeof(REAL_T), N * N);
    REAL_T *C_ref = calloc(sizeof(REAL_T), N * N);

    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            LOG_INFO("transa = %s, transb = %s\n", transa[i], transb[j]);

            LOG_INFO("alpha = %e\n", alpha);
            LOG_INFO("beta  = %e\n", beta);

            TYPED_FUNC(ref_random_matrix) (N, A);
            if (transa[i][0] == 'T')
            {
                TYPED_FUNC(transpose) (N, A, A_input);
            }
            else
            {
                TYPED_FUNC(copy_matrix) (N, A, A_input);
            }
            TYPED_FUNC(ref_random_matrix) (N, B);
            if (transb[j][0] == 'T')
            {
                TYPED_FUNC(transpose) (N, B, B_input);
            }
            else
            {
                TYPED_FUNC(copy_matrix) (N, B, B_input);
            }
            TYPED_FUNC(ref_random_matrix) (N, C);
            TYPED_FUNC(copy_matrix) (N, C, C_ref);

            LOG_INFO("A:\n");
            bml_print_dense_matrix(N, matrix_precision, dense_column_major,
                                   A_input, 0, N, 0, N);
            LOG_INFO("B:\n");
            bml_print_dense_matrix(N, matrix_precision, dense_column_major,
                                   B_input, 0, N, 0, N);
            LOG_INFO("C:\n");
            bml_print_dense_matrix(N, matrix_precision, dense_column_major, C,
                                   0, N, 0, N);

            TYPED_FUNC(ref_multiply) (N, A_input, B_input, C_ref, alpha, beta,
                                      0.0);
            LOG_INFO("C_ref:\n");
            bml_print_dense_matrix(N, matrix_precision, dense_column_major,
                                   C_ref, 0, N, 0, N);

            TYPED_FUNC(bml_gemm) (transa[i], transb[j], &N, &N, &N, &alpha, A,
                                  &N, B, &N, &beta, C, &N);
            LOG_INFO("C:\n");
            bml_print_dense_matrix(N, matrix_precision, dense_column_major, C,
                                   0, N, 0, N);

            if (TYPED_FUNC(compare_matrix) (N, matrix_precision, C, C_ref) !=
                0)
            {
                LOG_ERROR("matrix product incorrect\n");
                return -1;
            }
        }
    }
    LOG_INFO("test_bml_gemm passed\n");

    free(A);
    free(A_input);
    free(B);
    free(B_input);
    free(C);
    free(C_ref);

    return 0;
}
