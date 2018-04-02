#include "bml_introspection_ellpack.h"
#include "bml_types_ellpack.h"
#include "bml_logger.h"
#include "../../macros.h"
#include "../bml_types.h"
#include "../bml_introspection.h"
#include "../bml_logger.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>


/** Return the matrix precision.
 *
 * \param A The matrix.
 * \return The matrix precision.
 */
bml_matrix_precision_t
bml_get_precision_ellpack(
    const bml_matrix_ellpack_t * A)
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
 * \return The matrix distribution mode.
 */
bml_distribution_mode_t
bml_get_distribution_mode_ellpack(
    const bml_matrix_ellpack_t * A)
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
bml_get_N_ellpack(
    const bml_matrix_ellpack_t * A)
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
bml_get_M_ellpack(
    const bml_matrix_ellpack_t * A)
{
    if (A != NULL)
    {
        return A->M;
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
bml_get_row_bandwidth_ellpack(
    const bml_matrix_ellpack_t * A,
    const int i)
{
    return A->nnz[i];
}

/** Return the bandwidth of a matrix.
 *
 * \param A The bml matrix.
 * \return The bandwidth of row i.
 */
int
bml_get_bandwidth_ellpack(
    const bml_matrix_ellpack_t * A)
{
    int max_bandwidth = 0;
    for (int i = 0; i < A->N; i++)
    {
        max_bandwidth =
            (A->nnz[i] > max_bandwidth ? A->nnz[i] : max_bandwidth);
    }
    return max_bandwidth;
}


/** Return the sparsity of a matrix.
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
bml_get_sparsity_ellpack(
    const bml_matrix_ellpack_t * A,
    const double threshold)
{
    switch (bml_get_precision_ellpack(A))
    {
        case single_real:
            return bml_get_sparsity_ellpack_single_real(A, threshold);
            break;
        case double_real:
            return bml_get_sparsity_ellpack_double_real(A, threshold);
            break;
        case single_complex:
            return bml_get_sparsity_ellpack_single_complex(A, threshold);
            break;
        case double_complex:
            return bml_get_sparsity_ellpack_double_complex(A, threshold);
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
