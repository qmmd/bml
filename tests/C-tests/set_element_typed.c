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
    test_set_element) (
    const int N,
    const bml_matrix_type_t matrix_type,
    const bml_matrix_precision_t matrix_precision,
    const int M)
{
    bml_matrix_t *A = NULL;
    REAL_T *A_dense = NULL;
    REAL_T *B_dense = NULL;
    REAL_T *val = NULL;

    LOG_DEBUG("rel. tolerance = %e\n", REL_TOL);

    A_dense = calloc(N * N, sizeof(REAL_T));

    // Random dense matrix
    for (int i = 0; i < N * N; i++)
    {
        A_dense[i] = (REAL_T) (rand() / (double) RAND_MAX);
    }

    // Allocate a bml matrix
    A = bml_zero_matrix(matrix_type, matrix_precision, N, M, sequential);

    // Set a new element
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            val = &A_dense[ROWMAJOR(i, j, N, N)];
            bml_set_element_new(A, i, j, val);
        }
    }

    B_dense = bml_export_to_dense(A, dense_row_major);

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            REAL_T expected = A_dense[ROWMAJOR(i, j, N, N)];
            REAL_T actual = B_dense[ROWMAJOR(i, j, N, N)];
            double rel_diff = ABS((expected - actual) / expected);
            if (rel_diff > REL_TOL)
            {
                LOG_ERROR
                    ("bml_set_element_new: matrices are not identical; expected[%d] = %e, B[%d] = %e\n",
                     i, expected, i, actual);
                return -1;
            }
        }
    }

    // Free intermediate dense matrix before next test.
    bml_free_memory(B_dense);

    // Set an element that is not new
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            val = &A_dense[ROWMAJOR(i, j, N, N)];
            bml_set_element(A, i, j, val);
        }
    }

    B_dense = bml_export_to_dense(A, dense_row_major);

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            REAL_T expected = A_dense[ROWMAJOR(i, j, N, N)];
            REAL_T actual = B_dense[ROWMAJOR(i, j, N, N)];
            double rel_diff = ABS((expected - actual) / expected);
            if (rel_diff > REL_TOL)
            {
                LOG_ERROR
                    ("bml_set_element: matrices are not identical; expected[%d] = %e, B[%d] = %e\n",
                     i, expected, i, actual);
                return -1;
            }
        }
    }

    bml_free_memory(A_dense);
    bml_free_memory(B_dense);
    bml_deallocate(&A);

    LOG_INFO("get element test passed\n");

    return 0;
}
