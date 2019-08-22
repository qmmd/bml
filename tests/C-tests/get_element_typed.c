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
    test_get_element) (
    const int N,
    const bml_matrix_type_t matrix_type,
    const bml_matrix_precision_t matrix_precision,
    const int M)
{
    bml_matrix_t *A = NULL;
    REAL_T *A_dense = NULL;

    LOG_DEBUG("rel. tolerance = %e\n", REL_TOL);

    A_dense = calloc(N * N, sizeof(REAL_T));

    for (int i = 0; i < N * N; i++)
    {
        A_dense[i] = (REAL_T) (rand() / (double) RAND_MAX);
    }

    A = bml_import_from_dense(matrix_type, matrix_precision, dense_row_major,
                              N, M, A_dense, 0.0, sequential);

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            REAL_T *Aij = bml_get(A, i, j);
            REAL_T expected = A_dense[ROWMAJOR(i, j, N, N)];
            double rel_diff = ABS((expected - *Aij) / expected);
            if (rel_diff > REL_TOL)
            {
                LOG_ERROR
                    ("matrices are not identical; expected[%d] = %e, B[%d] = %e\n",
                     i, expected, i, *Aij);
                return -1;
            }
        }
    }
    bml_free_memory(A_dense);
    bml_deallocate(&A);

    LOG_INFO("get element test passed\n");

    return 0;
}
