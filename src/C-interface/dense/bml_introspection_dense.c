#include "../bml_logger.h"
#include "bml_introspection_dense.h"

#include <stdlib.h>

/** Return the matrix precision.
 *
 * \param A The matrix.
 * \return The matrix precision.
 */
bml_matrix_precision_t
bml_get_precision_dense(
    bml_matrix_dense_t * A)
{
    if (A != NULL)
    {
        return A->matrix_precision;
    }
    else
    {
        return precision_uninitialized;
    }
}

/** Return the matrix distribution mode.
 *
 * \param A The matrix.
 * \return The distribution mode.
 */
bml_distribution_mode_t
bml_get_distribution_mode_dense(
    bml_matrix_dense_t * A)
{
    if (A != NULL)
    {
        return A->distribution_mode;
    }
    else
    {
        return -1;
    }
}

/** Return the matrix size.
 *
 * \param A The matrix.
 * \return The matrix size.
 */
int
bml_get_N_dense(
    bml_matrix_dense_t * A)
{
    if (A != NULL)
    {
        return A->N;
    }
    else
    {
        return -1;
    }
}

/** Return the matrix parameter M.
 *
 * \param A The matrix.
 * \return The matrix parameter M.
 */
int
bml_get_M_dense(
    bml_matrix_dense_t * A)
{
    if (A != NULL)
    {
        return A->N;
    }
    else
    {
        return -1;
    }
}

/** Return the bandwidth of a row in the matrix.
 *
 * \param A The bml matrix.
 * \param i The row index.
 * \return The bandwidth of row i.
 */
int
bml_get_row_bandwidth_dense(
    bml_matrix_dense_t * A,
    int i)
{
    switch (bml_get_precision_dense(A))
    {
        case single_real:
            return bml_get_row_bandwidth_dense_single_real(A, i);
            break;
        case double_real:
            return bml_get_row_bandwidth_dense_double_real(A, i);
            break;
        case single_complex:
            return bml_get_row_bandwidth_dense_single_complex(A, i);
            break;
        case double_complex:
            return bml_get_row_bandwidth_dense_double_complex(A, i);
            break;
        case precision_uninitialized:
            LOG_ERROR("precision not initialized");
            break;
        default:
            LOG_ERROR("fatal logic error\n");
            break;
    }
    return -1;
}

/** Return the bandwidth of a matrix.
 *
 * \param A The bml matrix.
 * \return The bandwidth of row i.
 */
int
bml_get_bandwidth_dense(
    bml_matrix_dense_t * A)
{
    switch (bml_get_precision_dense(A))
    {
        case single_real:
            return bml_get_bandwidth_dense_single_real(A);
            break;
        case double_real:
            return bml_get_bandwidth_dense_double_real(A);
            break;
        case single_complex:
            return bml_get_bandwidth_dense_single_complex(A);
            break;
        case double_complex:
            return bml_get_bandwidth_dense_double_complex(A);
            break;
        case precision_uninitialized:
            LOG_ERROR("precision not initialized");
            break;
        default:
            LOG_ERROR("fatal logic error\n");
            break;
    }
    return -1;
}

/** Return the sparsity of a bml matrix.
 *
 *  Note that the the sparsity of a matrix is defined
 *  as NumberOfZeroes/N*N where N is the matrix dimension.
 *  The density of matrix A will be defined as 1-sparsity(A)
 *
 * \ingroup introspection_group_C
 *
 * \param A The bml matrix.
 * \param threshold The threshold used to compute the sparsity.
 * \return The sparsity of A.
 */
double
bml_get_sparsity_dense(
    bml_matrix_dense_t * A,
    double threshold)
{
    switch (bml_get_precision_dense(A))
    {
        case single_real:
            return bml_get_sparsity_dense_single_real(A, threshold);
            break;
        case double_real:
            return bml_get_sparsity_dense_double_real(A, threshold);
            break;
        case single_complex:
            return bml_get_sparsity_dense_single_complex(A, threshold);
            break;
        case double_complex:
            return bml_get_sparsity_dense_double_complex(A, threshold);
            break;
        case precision_uninitialized:
            LOG_ERROR("precision not initialized");
            break;
        default:
            LOG_ERROR("fatal logic error\n");
            break;
    }
    return -1;
}
