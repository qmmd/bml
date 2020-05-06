#include "bml_introspection.h"
#include "bml_logger.h"
#include "bml_types.h"
#include "dense/bml_introspection_dense.h"
#include "ellpack/bml_introspection_ellpack.h"
#include "ellsort/bml_introspection_ellsort.h"
#include "ellblock/bml_introspection_ellblock.h"
#include "csr/bml_introspection_csr.h"

#include <stdlib.h>

/** Returns the matrix type.
 *
 * If the matrix is not initialized yet, a type of "unitialized" is returned.
 *
 * \param A The matrix.
 * \return The matrix type.
 */
bml_matrix_type_t
bml_get_type(
    bml_matrix_t * A)
{
    bml_matrix_type_t *matrix_type = A;
    if (A != NULL)
    {
        return *matrix_type;
    }
    else
    {
        return type_uninitialized;
    }

}

/** Return the matrix precision.
 *
 * \param A The matrix.
 * \return The matrix precision.
 */
bml_matrix_precision_t
bml_get_precision(
    bml_matrix_t * A)
{
    switch (bml_get_type(A))
    {
        case type_uninitialized:
            return precision_uninitialized;
            break;
        case dense:
            return bml_get_precision_dense(A);
            break;
        case ellpack:
            return bml_get_precision_ellpack(A);
            break;
        case ellsort:
            return bml_get_precision_ellsort(A);
            break;
        case ellblock:
            return bml_get_precision_ellblock(A);
            break;
        case csr:
            return bml_get_precision_csr(A);
            break;
        default:
            LOG_ERROR("unknown precision");
            break;
    }
    return precision_uninitialized;
}

/** Return the matrix size.
 *
 * \param A The matrix.
 * \return The matrix size.
 */
int
bml_get_N(
    bml_matrix_t * A)
{
    switch (bml_get_type(A))
    {
        case type_uninitialized:
            return -1;
            break;
        case dense:
            return bml_get_N_dense(A);
            break;
        case ellpack:
            return bml_get_N_ellpack(A);
            break;
        case ellsort:
            return bml_get_N_ellsort(A);
            break;
        case ellblock:
            return bml_get_N_ellblock(A);
            break;
        case csr:
            return bml_get_N_csr(A);
            break;
        default:
            LOG_ERROR("unknown matrix type\n");
            break;
    }
    return -1;
}

/** Return the matrix parameter M.
 *
 * \param A The matrix.
 * \return The matrix parameter M.
 */
int
bml_get_M(
    bml_matrix_t * A)
{
    switch (bml_get_type(A))
    {
        case type_uninitialized:
            return 0;
            break;
        case dense:
            return bml_get_M_dense(A);
            break;
        case ellpack:
            return bml_get_M_ellpack(A);
            break;
        case ellsort:
            return bml_get_M_ellsort(A);
            break;
        case ellblock:
            return bml_get_M_ellblock(A);
            break;
        case csr:
            return bml_get_M_csr(A);
            break;
        default:
            LOG_ERROR("unknown matrix type\n");
            break;
    }
    return -1;
}

int
bml_get_NB(
    bml_matrix_t * A)
{
    switch (bml_get_type(A))
    {
        case type_uninitialized:
            return 0;
            break;
        case dense:
            return 1;
            break;
        case ellpack:
            return 1;
            break;
        case ellsort:
            return 1;
            break;
        case ellblock:
            return bml_get_NB_ellblock(A);
            break;
        case csr:
            return 1;
            break;
        default:
            LOG_ERROR("unknown matrix type\n");
            break;
    }
    return -1;
}

/** Return the bandwidth of a row in the matrix.
 *
 * \param A The bml matrix.
 * \param i The row index.
 * \return The bandwidth of row i.
 */
int
bml_get_row_bandwidth(
    bml_matrix_t * A,
    int i)
{
    if (i < 0 || i >= bml_get_N(A))
    {
        LOG_ERROR("row index %d outside of matrix\n", i);
    }
    switch (bml_get_type(A))
    {
        case dense:
            return bml_get_row_bandwidth_dense(A, i);
            break;
        case ellpack:
            return bml_get_row_bandwidth_ellpack(A, i);
            break;
        case ellsort:
            return bml_get_row_bandwidth_ellsort(A, i);
            break;
        case csr:
            return bml_get_row_bandwidth_csr(A, i);
            break;
        default:
            LOG_ERROR("unknown matrix type\n");
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
bml_get_bandwidth(
    bml_matrix_t * A)
{
    switch (bml_get_type(A))
    {
        case dense:
            return bml_get_bandwidth_dense(A);
            break;
        case ellpack:
            return bml_get_bandwidth_ellpack(A);
            break;
        case ellsort:
            return bml_get_bandwidth_ellsort(A);
            break;
        case csr:
            return bml_get_bandwidth_csr(A);
            break;
        default:
            LOG_ERROR("unknown matrix type\n");
            break;
    }
    return -1;
}

/** Return the distribution mode of a matrix.
 *
 * \param A The bml matrix.
 * \return The distibution mode of matrix A.
 */
bml_distribution_mode_t
bml_get_distribution_mode(
    bml_matrix_t * A)
{
    switch (bml_get_type(A))
    {
        case dense:
            return bml_get_distribution_mode_dense(A);
            break;
        case ellpack:
            return bml_get_distribution_mode_ellpack(A);
            break;
        case ellsort:
            return bml_get_distribution_mode_ellsort(A);
            break;
        case csr:
            return bml_get_distribution_mode_csr(A);
            break;
        default:
            LOG_ERROR("unknown matrix type in bml_get_distribution_mode\n");
            break;
    }
    return -1;
}

/** Return the sparsity of a matrix.
 *
 * \param A The bml matrix.
 * \param threshold The threshold used to compute the sparsity.
 * \return The sparsity of matrix A.
 */
double
bml_get_sparsity(
    bml_matrix_t * A,
    double threshold)
{
    switch (bml_get_type(A))
    {
        case dense:
            return bml_get_sparsity_dense(A, threshold);
            break;
        case ellpack:
            return bml_get_sparsity_ellpack(A, threshold);
            break;
        case ellsort:
            return bml_get_sparsity_ellsort(A, threshold);
            break;
        case ellblock:
            return bml_get_sparsity_ellblock(A, threshold);
            break;
        case csr:
            return bml_get_sparsity_csr(A, threshold);
            break;
        default:
            LOG_ERROR("unknown matrix type in bml_get_sparsity\n");
            break;
    }
    return -1;
}
