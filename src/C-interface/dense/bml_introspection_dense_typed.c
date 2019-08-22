#include "../../macros.h"
#include "../../typed.h"
#include "../bml_export.h"
#include "bml_export_dense.h"
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
    bml_matrix_dense_t * A,
    int i)
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
    bml_matrix_dense_t * A)
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
    bml_matrix_dense_t * A,
    double threshold)
{

#ifdef BML_USE_MAGMA
    REAL_T *A_matrix = bml_export_to_dense(A, dense_row_major);
#else
    REAL_T *A_matrix = A->matrix;
#endif
    int nnzs = 0;
    int N = A->N;
    double sparsity;

    nnzs = 0;
    for (int i = 0; i < A->N; i++)
    {
        for (int j = 0; j < A->N; j++)
        {
            if (is_above_threshold
                (A_matrix[ROWMAJOR(i, j, A->N, A->N)], threshold))
            {
                nnzs++;
            }
        }
    }
#ifdef BML_USE_MAGMA
    free(A_matrix);
#endif
    sparsity = (1.0 - (double) nnzs / ((double) (N * N)));

    return sparsity;
}
