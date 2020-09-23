#include "../../macros.h"
#include "../bml_introspection.h"
#include "../bml_logger.h"
#include "../bml_logger.h"
#include "../bml_types.h"
#include "bml_introspection_csr.h"
#include "bml_types_csr.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>


/** Return the matrix precision.
 *
 * \param A The matrix.
 * \return The matrix precision.
 */
bml_matrix_precision_t
bml_get_precision_csr(
    bml_matrix_csr_t * A)
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
bml_get_distribution_mode_csr(
    bml_matrix_csr_t * A)
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
bml_get_N_csr(
    bml_matrix_csr_t * A)
{
    if (A != NULL)
    {
        return A->N_;
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
bml_get_M_csr(
    bml_matrix_csr_t * A)
{
    if (A != NULL)
    {
        int max_alloc = 0;
        const int N = A->N_;
        for (int i = 0; i < N; i++)
        {
            max_alloc =
                (csr_row_alloc_size(A->data_[i]) >
                 max_alloc ? csr_row_alloc_size(A->data_[i]) : max_alloc);
        }

        return max_alloc;
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
bml_get_row_bandwidth_csr(
    bml_matrix_csr_t * A,
    int i)
{
    return csr_row_NNZ(A->data_[i]);
}

/** Return the bandwidth of a matrix.
 *
 * \param A The bml matrix.
 * \return The bandwidth of A.
 */
int
bml_get_bandwidth_csr(
    bml_matrix_csr_t * A)
{
    int max_bandwidth = 0;
    const int N = A->N_;
#pragma omp parallel for reduction (max:max_bandwidth)
    for (int i = 0; i < N; i++)
    {
        max_bandwidth =
            (csr_row_NNZ(A->data_[i]) >
             max_bandwidth ? csr_row_NNZ(A->data_[i]) : max_bandwidth);
    }
    return max_bandwidth;
}


/** Return the sparsity of a matrix.
 *
 *  Note that the the sparsity of a matrix is defined
 *  as (1 - NumberOfZeroes/N*N) where N is the matrix dimension.
 *  The density of matrix A will be defined as 1-sparsity(A)
 *
 * \ingroup introspection_group_C
 *
 * \param A The bml matrix.
 * \param threshold The threshold used to compute the sparsity.
 * \return The sparsity of A.
 */
double
bml_get_sparsity_csr(
    bml_matrix_csr_t * A,
    double threshold)
{
    switch (bml_get_precision_csr(A))
    {
        case single_real:
            return bml_get_sparsity_csr_single_real(A, threshold);
            break;
        case double_real:
            return bml_get_sparsity_csr_double_real(A, threshold);
            break;
        case single_complex:
            return bml_get_sparsity_csr_single_complex(A, threshold);
            break;
        case double_complex:
            return bml_get_sparsity_csr_double_complex(A, threshold);
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
