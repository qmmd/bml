#include "bml.h"
#include "../typed.h"
#include "../macros.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#if defined(SINGLE_REAL) || defined(SINGLE_COMPLEX)
#define REL_TOL 1e-6
#else
#define REL_TOL 1e-12
#endif

int TYPED_FUNC(
    test_bml_add) (
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

    int result = 0;

    double norm;
    double expected_norm;

    double alpha = 1.2;
    double beta = 0.8;
    double threshold = 0.0;

    LOG_DEBUG("rel. tolerance = %e\n", REL_TOL);

    A = bml_random_matrix(matrix_type, matrix_precision, N, M, sequential);
    A_dense = bml_export_to_dense(A, dense_row_major);

    LOG_INFO("A_dense:\n");
    bml_print_dense_matrix(N, matrix_precision, dense_row_major, A_dense, 0,
                           N, 0, N);
    C = bml_random_matrix(matrix_type, matrix_precision, N, M, sequential);
    C_dense = bml_export_to_dense(C, dense_row_major);

    LOG_INFO("C_dense:\n");
    bml_print_dense_matrix(N, matrix_precision, dense_row_major, C_dense, 0,
                           N, 0, N);
    LOG_INFO("Testing bml_add()\n");

    B = bml_copy_new(A);
    bml_add(B, C, alpha, beta, threshold);
    B_dense = bml_export_to_dense(B, dense_row_major);

    LOG_INFO("B_dense:\n");
    bml_print_dense_matrix(N, matrix_precision, dense_row_major, B_dense, 0,
                           N, 0, N);
    for (int i = 0; i < N * N; i++)
    {
        double expected = alpha * A_dense[i] + beta * C_dense[i];
        double rel_diff = ABS((expected - B_dense[i]) / expected);
        if (rel_diff > REL_TOL)
        {
            LOG_INFO
                ("matrices are not identical; expected[%d] = %e, B[%d] = %e\n",
                 i, expected, i, B_dense[i]);
            result = -1;
        }
    }
    bml_free_memory(B_dense);
    bml_deallocate(&B);

    LOG_INFO("Testing bml_add_norm()\n");

    B = bml_copy_new(A);
    norm = bml_add_norm(B, C, alpha, beta, threshold);
    B_dense = bml_export_to_dense(B, dense_row_major);

    LOG_INFO("B_dense:\n");
    bml_print_dense_matrix(N, matrix_precision, dense_row_major, B_dense, 0,
                           N, 0, N);
    expected_norm = 0;
    for (int i = 0; i < N * N; i++)
    {
        expected_norm += C_dense[i] * C_dense[i];
        double expected = alpha * A_dense[i] + beta * C_dense[i];
        double rel_diff = ABS((expected - B_dense[i]) / expected);
        if (rel_diff > REL_TOL)
        {
            LOG_INFO
                ("matrices are not identical; expected[%d] = %e, B[%d] = %e\n",
                 i, expected, i, B_dense[i]);
            result = -1;
        }
    }
    if (ABS(expected_norm - norm) / expected_norm > REL_TOL)
    {
        LOG_INFO("norm mismatch: expected = %e, norm = %e\n", expected_norm, norm);
        result = -1;
    }
    bml_free_memory(B_dense);
    bml_deallocate(&B);

    LOG_INFO("Testing bml_add_identity()\n");

    B = bml_copy_new(A);
    bml_add_identity(B, beta, threshold);
    B_dense = bml_export_to_dense(B, dense_row_major);

    LOG_INFO("B_dense:\n");
    bml_print_dense_matrix(N, matrix_precision, dense_row_major, B_dense, 0,
                           N, 0, N);
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            double expected = 0;
            if (i == j)
            {
                expected = A_dense[ROWMAJOR(i, i, N, N)] + beta;
            }
            else
            {
                expected = A_dense[ROWMAJOR(i, j, N, N)];
            }

            double rel_diff =
                ABS((expected - B_dense[ROWMAJOR(i, j, N, N)]) / expected);
            if (rel_diff > REL_TOL)
            {
                LOG_INFO
                    ("matrices are not identical; expected[%d] = %e, B[%d] = %e\n",
                     i, expected, i, B_dense[ROWMAJOR(i, j, N, N)]);
                result = -1;
            }
        }
    }
    bml_free_memory(B_dense);
    bml_deallocate(&B);

    LOG_INFO("Testing bml_scale_add_identity()\n");

    B = bml_copy_new(A);
    bml_scale_add_identity(B, alpha, beta, threshold);
    B_dense = bml_export_to_dense(B, dense_row_major);

    LOG_INFO("B_dense:\n");
    bml_print_dense_matrix(N, matrix_precision, dense_row_major, B_dense, 0,
                           N, 0, N);
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            double expected = 0;
            if (i == j)
            {
                expected = alpha * A_dense[ROWMAJOR(i, i, N, N)] + beta;
            }
            else
            {
                expected = alpha * A_dense[ROWMAJOR(i, j, N, N)];
            }

            double rel_diff =
                ABS((expected - B_dense[ROWMAJOR(i, j, N, N)]) / expected);
            if (rel_diff > REL_TOL)
            {
                LOG_INFO
                    ("matrices are not identical; expected[%d] = %e, B[%d] = %e\n",
                     i, expected, i, B_dense[ROWMAJOR(i, j, N, N)]);
                result = -1;
            }
        }
    }
    bml_free_memory(B_dense);
    bml_deallocate(&B);

    bml_free_memory(A_dense);
    bml_free_memory(C_dense);

    bml_deallocate(&A);
    bml_deallocate(&C);

    if (result == 0)
    {
        LOG_INFO("add matrix test passed\n");
    }
    else
    {
        LOG_INFO("add matrix test failed\n");
    }

    return result;
}
