#include "../../macros.h"
#include "bml_getters_dense.h"
#include "../bml_introspection.h"
#include "bml_types_dense.h"
#include "../typed.h"
#include "../bml_logger.h"

#include <complex.h>
#include <stdlib.h>

REAL_T *TYPED_FUNC(
    bml_get_dense) (
    const bml_matrix_dense_t * A,
    const int i,
    const int j)
{
    int N = bml_get_N(A);
    REAL_T *A_matrix = A->matrix;

    return &A_matrix[ROWMAJOR(i, j, N, N)];
}

void *TYPED_FUNC(
    bml_get_row_dense) (
    bml_matrix_dense_t * A,
    const int i)
{
    int N = bml_get_N(A);
    REAL_T *A_matrix = A->matrix;
    REAL_T *row = calloc(N, sizeof(REAL_T));

    LOG_DEBUG("here\n");

    for (int j = 0; j < N; j++)
    {
        row[j] = A_matrix[ROWMAJOR(i, j, N, N)];
    }

    return row;
}

void *TYPED_FUNC(
    bml_get_diagonal_dense) (
    bml_matrix_dense_t * A)
{
    int N = bml_get_N(A);
    REAL_T *A_matrix = A->matrix;
    REAL_T *diagonal = calloc(N, sizeof(REAL_T));

    for (int j = 0; j < N; j++)
    {
        diagonal[j] = A_matrix[ROWMAJOR(j, j, N, N)];
    }

    return diagonal;
}
