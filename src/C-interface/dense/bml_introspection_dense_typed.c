#include "../../macros.h"
#include "../../typed.h"
#include "bml_introspection_dense.h"

#include <assert.h>
#include <complex.h>
#include <math.h>
#include <stdlib.h>

/** Return the bandwidth of a row in the matrix.
 *
 * \param A The bml matrix.
 * \param i The row index.
 * \return The bandwidth of row i.
 */
int TYPED_FUNC(
    bml_get_row_bandwidth_dense) (
    const bml_matrix_dense_t * A,
    const int i)
{
    assert(A != NULL);

    REAL_T *A_matrix = A->matrix;
    int bandwidth = 0;

    for (int j = 0; j < A->N; j++)
    {
        if (is_above_threshold(A_matrix[ROWMAJOR(i, j, A->N, A->N)], 0))
        {
            bandwidth++;
        }
    }
    return bandwidth;
}

/** Return the bandwidth of a matrix.
 *
 * \param A The bml matrix.
 * \return The bandwidth of the matrix.
 */
int TYPED_FUNC(
    bml_get_bandwidth_dense) (
    const bml_matrix_dense_t * A)
{
    assert(A != NULL);

    REAL_T *A_matrix = A->matrix;
    int bandwidth = 0;
    int row_bandwidth = 0;

    for (int i = 0; i < A->N; i++)
    {
        row_bandwidth = 0;
        for (int j = 0; j < A->N; j++)
        {
            if (is_above_threshold(A_matrix[ROWMAJOR(i, j, A->N, A->N)], 0))
            {
                row_bandwidth++;
            }
        }
        if (row_bandwidth > bandwidth)
        {
            bandwidth = row_bandwidth;
        }
    }
    return bandwidth;
}

/** Return the sparsity of a matrix.
 *
 * \param A The bml matrix.
 * \param threshold The threshold used to compute the sparsity.
 * \return The sparsity of a matrix.
 */
double TYPED_FUNC(
    bml_get_sparsity_dense) (
    const bml_matrix_dense_t * A,
    const double threshold)
{

    REAL_T *A_matrix = A->matrix;
    int nnzs = 0;
    int N = A->N;
    double sparsity;

    nnzs = 0;
    for (int i = 0; i < A->N; i++)
    {
        for (int j = 0; j < A->N; j++)
        {
            if (ABS(A_matrix[ROWMAJOR(i, j, A->N, A->N)]) > threshold)
            {
                nnzs++;
            }
        }
    }

    sparsity = (1.0 - (double) nnzs / ((double) (N * N)));

    return sparsity;
}
