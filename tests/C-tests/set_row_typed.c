#include "bml.h"
#include "../typed.h"
#include "bml_getters.h"
#include "bml_setters.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

int TYPED_FUNC(
    test_set_row) (
    const int N,
    const bml_matrix_type_t matrix_type,
    const bml_matrix_precision_t matrix_precision,
    const int M)
{
    bml_matrix_t *A = NULL;
    REAL_T *A_row = NULL;
    REAL_T *B_row = NULL;

    // Create a row
    switch (matrix_precision)
    {
        case single_real:
            A_row = calloc(N, sizeof(float));
            break;
        case double_real:
            A_row = calloc(N, sizeof(double));
            break;
        case single_complex:
            A_row = calloc(N, sizeof(float complex));
            break;
        case double_complex:
            A_row = calloc(N, sizeof(double complex));
            break;
        default:
            LOG_DEBUG("matrix_precision is not set");
            break;
    }

    for (int i = 0; i < N; i++)
    {
        A_row[i] = i;
    }

    A = bml_zero_matrix(matrix_type, matrix_precision, N, M, sequential);

    // Set the second row
    bml_set_row(A, 2, A_row, 1.0e-10);

    // Retrive the second row
    B_row = bml_get_row(A, 2);

    for (int i = 0; i < N; i++)
    {
        if (ABS(A_row[i] - B_row[i]) > 1e-12)
        {
            LOG_ERROR("bml_set_row and/or bml_get_row are corrupted\n");
            return -1;
        }
    }

    bml_print_bml_matrix(A, 0, N, 0, N);

    bml_deallocate(&A);
    free(A_row);
    free(B_row);

    LOG_INFO("bml_set_row passed\n");

    return 0;
}
