#include "bml.h"
#include "../macros.h"
#include "../typed.h"
#include "ellblock/bml_allocate_ellblock.h"
#include "bml_utilities.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

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
    REAL_T * A,
    REAL_T * B)
{
    LOG_INFO("Compare matrices\n");
    for (int i = 0; i < N * N; i++)
    {
        if (ABS(A[i] - B[i]) > ABS_TOL)
        {
            LOG_INFO("First matrix\n");
            bml_print_dense_matrix(N, matrix_precision, dense_row_major, A, 0,
                                   N, 0, N);
            LOG_INFO("Second matrix\n");
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

void TYPED_FUNC(
    setup_bsizes) (
    const int N,
    const int M)
{
    int bsizes[N];
    int count = 0;
    int nb = 0;
    for (int i = 0; i < N; i++)
    {
        int bsize = (3 + i) % 6 + 1;
        count += bsize;
        if (count > N)
        {
            bsize -= (count - N);
            count = N;
        }
        bsizes[i] = bsize;
        if (count == N)
        {
            nb = i + 1;
            break;
        }
    }
    bml_set_block_sizes(bsizes, nb, nb * M / N);
    for (int i = 0; i < nb; i++)
    {
        printf("bsizes[%d] = %d\n", i, bsizes[i]);
    }

}

int TYPED_FUNC(
    test_multiply) (
    const int N,
    const bml_matrix_type_t matrix_type,
    const bml_matrix_precision_t matrix_precision,
    const int M)
{
    // set block sizes for ellblock tests
    // (unused by other tests)
#ifdef DO_MPI
    if (bml_getNRanks() > 1)
    {
        int p2 = bml_getNRanks();
        int p = bml_sqrtint(p2);
        TYPED_FUNC(setup_bsizes) (N / p, M / p);
    }
    else
#endif
    if (matrix_type == ellblock)
    {
        LOG_INFO("ellblock type\n");
        TYPED_FUNC(setup_bsizes) (N, M);
    }

    bml_distribution_mode_t distrib_mode = sequential;
#ifdef DO_MPI
    if (bml_getNRanks() > 1)
    {
        LOG_INFO("Use distributed matrix\n");
        distrib_mode = distributed;
    }
#endif

    const double alpha = 1.2;
    const double beta = 0.8;
    const double threshold = 0.0;

    bml_matrix_t *A =
        bml_random_matrix(matrix_type, matrix_precision, N, M, distrib_mode);
    bml_matrix_t *B =
        bml_random_matrix(matrix_type, matrix_precision, N, M, distrib_mode);
    bml_matrix_t *C =
        bml_random_matrix(matrix_type, matrix_precision, N, M, distrib_mode);

    REAL_T *A_dense = bml_export_to_dense(A, dense_row_major);
    REAL_T *B_dense = bml_export_to_dense(B, dense_row_major);
    REAL_T *C_dense = bml_export_to_dense(C, dense_row_major);
    REAL_T *D_dense = bml_export_to_dense(C, dense_row_major);

    LOG_INFO("bml_multiply\n");
    bml_multiply(A, B, C, alpha, beta, threshold);

    LOG_INFO("bml_export_to_dense\n");
    REAL_T *E_dense = bml_export_to_dense(C, dense_row_major);

    if (bml_getMyRank() == 0)
    {
        TYPED_FUNC(ref_multiply) (N, A_dense, B_dense, D_dense, alpha, beta,
                                  threshold);
        bml_print_dense_matrix(N, matrix_precision, dense_row_major, E_dense,
                               0, N, 0, N);
        bml_print_dense_matrix(N, matrix_precision, dense_row_major, D_dense,
                               0, N, 0, N);

        if (TYPED_FUNC(compare_matrix) (N, matrix_precision, D_dense, E_dense)
            != 0)
        {
            LOG_ERROR("matrix product incorrect\n");
            return -1;
        }
        LOG_INFO("multiply matrix test passed\n");
    }

    bml_deallocate(&A);
    bml_deallocate(&B);
    bml_deallocate(&C);

    if (bml_getMyRank() == 0)
    {
        bml_free_memory(A_dense);
        bml_free_memory(B_dense);
        bml_free_memory(C_dense);
        bml_free_memory(D_dense);
        bml_free_memory(E_dense);
    }
    return 0;
}
