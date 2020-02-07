#include "bml.h"
#include "../typed.h"

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
    test_trace_mult) (
    const int N,
    const bml_matrix_type_t matrix_type,
    const bml_matrix_precision_t matrix_precision,
    const int M)
{
    bml_matrix_t *A = NULL;
    bml_matrix_t *B = NULL;

    double trace, trace_mult, diff, threshold;
    threshold = 0.0;

    A = bml_random_matrix(matrix_type, matrix_precision, N, M, sequential);
    B = bml_random_matrix(matrix_type, matrix_precision, N, M, sequential);

    bml_matrix_t *A_t = bml_transpose_new(A);
    bml_matrix_t *B_t = bml_transpose_new(B);

    bml_add(A, A_t, 1, 1, 0);
    bml_add(B, B_t, 1, 1, 0);

    printf("A\n");
    bml_print_bml_matrix(A, 0, N, 0, N);

    printf("B\n");
    bml_print_bml_matrix(B, 0, N, 0, N);

    trace = 0;
    for (int i = 0; i < N; i++)
    {
        for (int k = 0; k < N; k++)
        {
            trace +=
                (*(REAL_T *) bml_get(A, i, k)) *
                (*(REAL_T *) bml_get(B, k, i));
        }
    }

    trace_mult = bml_trace_mult(A, B);

    diff = fabs(trace - trace_mult) / (double) (N * N);
    if (diff > REL_TOL)
    {
        LOG_ERROR
            ("trace_mult is not correct; trace = %e and not %e, diff / N^2 = %e\n",
             trace_mult, trace, diff);
        return -1;
    }
    LOG_INFO("test_trace_mult passed (diff / N^2 = %e)\n", diff);

    bml_deallocate(&A);
    bml_deallocate(&A_t);
    bml_deallocate(&B);
    bml_deallocate(&B_t);

    return 0;
}
