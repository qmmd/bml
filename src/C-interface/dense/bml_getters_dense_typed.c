#ifdef BML_USE_MAGMA
#include "magma_v2.h"
#endif

#include "../../macros.h"
#include "../../typed.h"
#include "../bml_export.h"
#include "../bml_introspection.h"
#include "../bml_logger.h"
#include "bml_export_dense.h"
#include "bml_getters_dense.h"
#include "bml_types_dense.h"

#include <complex.h>
#include <stdlib.h>

void *TYPED_FUNC(
    bml_get_dense) (
    bml_matrix_dense_t * A,
    int i,
    int j)
{
    int N = bml_get_N(A);

    if (N < 0)
    {
        LOG_ERROR("A is not initialized\n");
    }

#ifdef BML_USE_MAGMA
    static REAL_T value = 0.;
    MAGMA_T tmp;
    MAGMA(getvector) (1, (MAGMA_T *) A->matrix + i * A->ld + j, 1,
                      &tmp, 1, A->queue);
    value = MAGMACOMPLEX(REAL) (tmp);
    return &value;
#else
    REAL_T *A_matrix = A->matrix;

    return &A_matrix[ROWMAJOR(i, j, N, N)];
#endif
}

void *TYPED_FUNC(
    bml_get_row_dense) (
    bml_matrix_dense_t * A,
    int i)
{
    int N = bml_get_N(A);

    if (N < 0)
    {
        LOG_ERROR("A is not initialized\n");
    }

#ifdef BML_USE_MAGMA
    MAGMA_T *A_matrix = A->matrix;
    MAGMA_T *row = calloc(N, sizeof(MAGMA_T));

    MAGMA(getvector) (N, (MAGMA_T *) A->matrix + i * A->ld, 1,
                      row, 1, A->queue);
#else
    REAL_T *A_matrix = A->matrix;
    REAL_T *row = calloc(N, sizeof(REAL_T));

    for (int j = 0; j < N; j++)
    {
        row[j] = A_matrix[ROWMAJOR(i, j, N, N)];
    }
#endif
    return row;
}

void *TYPED_FUNC(
    bml_get_diagonal_dense) (
    bml_matrix_dense_t * A)
{
    int N = bml_get_N(A);

    if (N < 0)
    {
        LOG_ERROR("A is not initialized\n");
    }

#ifdef BML_USE_MAGMA
    REAL_T *A_matrix = bml_export_to_dense(A, dense_row_major);
#else
    REAL_T *A_matrix = A->matrix;
#endif
    REAL_T *diagonal = calloc(N, sizeof(REAL_T));

    for (int j = 0; j < N; j++)
    {
        diagonal[j] = A_matrix[ROWMAJOR(j, j, N, N)];
    }
#ifdef BML_USE_MAGMA
    free(A_matrix);
#endif

    return diagonal;
}
