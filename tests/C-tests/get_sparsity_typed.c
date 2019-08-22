#include "bml.h"
#include "../macros.h"
#include "../typed.h"
#include "bml_introspection.h"
#include "bml_export.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

int TYPED_FUNC(
    test_get_sparsity) (
    const int N,
    const bml_matrix_type_t matrix_type,
    const bml_matrix_precision_t matrix_precision,
    const int M)
{
    bml_matrix_t *A = NULL;
    double sparsity;
    double sparsity_ref;
    int nnzs = 0;
    double threshold = 0.1;
    REAL_T *A_dense;

    A = bml_random_matrix(matrix_type, matrix_precision, N, M, sequential);

    sparsity = bml_get_sparsity(A, threshold);

    // Form bml to dense_array
    A_dense = bml_export_to_dense(A, dense_row_major);

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if (ABS(A_dense[ROWMAJOR(i, j, N, N)]) > threshold)
            {
                nnzs++;
            }
        }
    }

    sparsity_ref = (1.0 - (double) nnzs / ((double) (N * N)));

    if (fabs(sparsity - sparsity_ref) > 1e-12)
    {
        LOG_ERROR("bml_get_sparsity is corrupted\n");
        return -1;
    }

    printf("Sparsity = %f\n", sparsity);

    bml_deallocate(&A);
    free(A_dense);

    LOG_INFO("bml_get_sparsity passed\n");

    return 0;
}
