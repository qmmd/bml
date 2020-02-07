#include "bml.h"
#include "../typed.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

int TYPED_FUNC(
    test_template) (
    const int N,
    const bml_matrix_type_t matrix_type,
    const bml_matrix_precision_t matrix_precision,
    const int M)
{
    bml_matrix_t *A = NULL;
    REAL_T *array = NULL;

    // Add code here ...
    A = bml_random_matrix(matrix_type, matrix_precision, N, M, sequential);

    // Don't forget to release the memory ...
    bml_deallocate(&A);
    free(array);

    LOG_INFO("test_template passed\n");

    return 0;
}
