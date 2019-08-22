#include "bml_getters_ellsort.h"
#include "../bml_introspection.h"
#include "../bml_logger.h"
#include "bml_types_ellsort.h"
#include "../../macros.h"
#include "../../typed.h"

#include <complex.h>
#include <stdlib.h>

/** Return a single matrix element.
 *
 * \param A The bml matrix
 * \param i The row index
 * \param j The column index
 * \return The matrix element
 */
void *TYPED_FUNC(
    bml_get_ellsort) (
    bml_matrix_ellsort_t * A,
    int i,
    int j)
{
    static REAL_T MINUS_ONE = -1;
    static REAL_T ZERO = 0;
    REAL_T *A_value = (REAL_T *) A->value;

    if (i < 0 || i >= A->N)
    {
        LOG_ERROR("row index out of bounds\n");
        return &MINUS_ONE;
    }
    if (j < 0 || j >= A->N)
    {
        LOG_ERROR("column index out of bounds\n");
        return &MINUS_ONE;
    }
    for (int j_index = 0; j_index < A->nnz[i]; j_index++)
    {
        if (A->index[ROWMAJOR(i, j_index, A->N, A->M)] == j)
        {
            return &A_value[ROWMAJOR(i, j_index, A->N, A->M)];
        }
    }
    return &ZERO;
}

/** Get row i of matrix A.
 *
 *  \ingroup getters
 *
 *  \param A The matrix which takes row i
 *  \param i The index of the row to get
 *  \param row Array to copy the row
 *
 */
void *TYPED_FUNC(
    bml_get_row_ellsort) (
    bml_matrix_ellsort_t * A,
    int i)
{
    int ll;
    int A_N = A->N;
    int A_M = A->M;
    REAL_T *A_value = (REAL_T *) A->value;
    int *A_index = A->index;
    int *A_nnz = A->nnz;
    REAL_T *row = calloc(A_N, sizeof(REAL_T));

    for (int i = 0; i < A_N; i++)
    {
        row[i] = 0.0;
    }

    for (int j = 0; j < A_nnz[i]; j++)
    {
        ll = A_index[ROWMAJOR(i, j, A_N, A_M)];
        if (ll >= 0)
        {
            row[ll] = A_value[ROWMAJOR(i, j, A_N, A_M)];
        }
    }

    return row;
}

/** Get the diagonal of matrix A.
 *
 *  \ingroup getters
 *
 *  \param A The matrix which takes row i
 *  \param Diagonal Array to copy the diagonal
 *
 */
void *TYPED_FUNC(
    bml_get_diagonal_ellsort) (
    bml_matrix_ellsort_t * A)
{
    int A_N = A->N;
    int A_M = A->M;
    REAL_T *A_value = (REAL_T *) A->value;
    int *A_index = A->index;
    int *A_nnz = A->nnz;
    REAL_T *diagonal = calloc(A_N, sizeof(REAL_T));

    for (int i = 0; i < A_N; i++)
    {
        diagonal[i] = 0.0;
        for (int j = 0; j < A_nnz[i]; j++)
        {
            if (A_index[ROWMAJOR(i, j, A_N, A_M)] == i)
            {
                diagonal[i] = A_value[ROWMAJOR(i, j, A_N, A_M)];
            }
        }
    }

    return diagonal;
}
