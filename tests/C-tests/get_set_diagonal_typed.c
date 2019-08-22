#include "bml.h"
#include "../typed.h"
#include "bml_getters.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

int TYPED_FUNC(
    test_get_set_diagonal) (
    const int N,
    const bml_matrix_type_t matrix_type,
    const bml_matrix_precision_t matrix_precision,
    const int M)
{
    bml_matrix_t *A = NULL;
    REAL_T *A_diagonal = NULL;
    REAL_T *B_diagonal = NULL;

    // Create a diagonal
    switch (matrix_precision)
    {
        case single_real:
            A_diagonal = calloc(N, sizeof(float));
            break;
        case double_real:
            A_diagonal = calloc(N, sizeof(double));
            break;
        case single_complex:
            A_diagonal = calloc(N, sizeof(float complex));
            break;
        case double_complex:
            A_diagonal = calloc(N, sizeof(double complex));
            break;
        default:
            LOG_DEBUG("matrix_precision is not set");
            break;
    }


    for (int i = 0; i < N; i++)
    {
        A_diagonal[i] = i;
    }

    A = bml_random_matrix(matrix_type, matrix_precision, N, M, sequential);

    bml_set_diagonal(A, A_diagonal, 0.01);

    bml_print_bml_matrix(A, 0, N, 0, N);

    B_diagonal = bml_get_diagonal(A);

    for (int i = 0; i < N; i++)
    {
        if (ABS(A_diagonal[i] - B_diagonal[i]) > 1e-12)
        {
            LOG_ERROR
                ("bml_get_diagonal and/or bml_set_diagonal are corrupted\n");
            return -1;
        }
    }

    free(A_diagonal);
    free(B_diagonal);
    bml_deallocate(&A);

    LOG_INFO("bml_get_diagonal and bml_set_diagonal test passed\n");

    return 0;
}
