#include "../../macros.h"
#include "../../typed.h"
#include "../bml_introspection.h"
#include "bml_setters_dense.h"
#include "bml_types_dense.h"
#include "../bml_utilities.h"

#include <complex.h>
#include <stdio.h>

void TYPED_FUNC(
    bml_set_element_dense) (
    bml_matrix_dense_t * A,
    const int i,
    const int j,
    const void *value)
{
    int N = bml_get_N(A);
    REAL_T *A_matrix = A->matrix;
    A_matrix[ROWMAJOR(i, j, N, N)] = *((REAL_T *) value);
}

void TYPED_FUNC(
    bml_set_row_dense) (
    bml_matrix_dense_t * A,
    const int i,
    const REAL_T * row)
{
    int N = bml_get_N(A);
    REAL_T *A_matrix = A->matrix;
    for (int j = 0; j < N; j++)
    {
        A_matrix[ROWMAJOR(i, j, N, N)] = row[j];
    }
}

void TYPED_FUNC(
    bml_set_diagonal_dense) (
    bml_matrix_dense_t * A,
    const REAL_T * diagonal)
{
    int N = bml_get_N(A);
    REAL_T *A_matrix = A->matrix;
    for (int j = 0; j < N; j++)
    {
        A_matrix[ROWMAJOR(j, j, N, N)] = diagonal[j];
    }
}
