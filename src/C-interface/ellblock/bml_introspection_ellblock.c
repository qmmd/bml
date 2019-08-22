#include "../../macros.h"
#include "../bml_introspection.h"
#include "../bml_logger.h"
#include "../bml_logger.h"
#include "../bml_types.h"
#include "bml_introspection_ellblock.h"
#include "bml_types_ellblock.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>


/** Return the matrix precision.
 *
 * \param A The matrix.
 * \return The matrix precision.
 */
bml_matrix_precision_t
bml_get_precision_ellblock(
    bml_matrix_ellblock_t * A)
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
bml_get_distribution_mode_ellblock(
    bml_matrix_ellblock_t * A)
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
bml_get_N_ellblock(
    bml_matrix_ellblock_t * A)
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
bml_get_M_ellblock(
    bml_matrix_ellblock_t * A)
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

int
bml_get_NB_ellblock(
    bml_matrix_ellblock_t * A)
{
    if (A != NULL)
    {
        return A->NB;
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
bml_get_row_bandwidth_ellblock(
    bml_matrix_ellblock_t * A,
    int i)
{
    //determine block index and index within block
    int ib = 0;
    int ii = i;
    while (ii >= A->bsize[ib])
    {
        ii -= A->bsize[ib];
        ib++;
    }
    int row_bandwidth = 0;
    for (int jp = 0; jp < A->nnzb[ib]; jp++)
    {
        int jb = A->indexb[ROWMAJOR(ib, jp, A->NB, A->MB)];
        row_bandwidth += A->bsize[jb];
    }
    return row_bandwidth;
}

/** Return the bandwidth of a matrix.
 *
 * \param A The bml matrix.
 * \return The bandwidth of row i.
 */
int
bml_get_bandwidth_ellblock(
    bml_matrix_ellblock_t * A)
{
    int max_bandwidth = 0;
    int offset = 0;
    for (int ib = 0; ib < A->NB; ib++)
    {
        int bw = bml_get_row_bandwidth_ellblock(A, offset);
        max_bandwidth = (bw > max_bandwidth ? bw : max_bandwidth);
        offset += A->bsize[ib];
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
bml_get_sparsity_ellblock(
    bml_matrix_ellblock_t * A,
    double threshold)
{
    switch (bml_get_precision_ellblock(A))
    {
        case single_real:
            return bml_get_sparsity_ellblock_single_real(A, threshold);
            break;
        case double_real:
            return bml_get_sparsity_ellblock_double_real(A, threshold);
            break;
        case single_complex:
            return bml_get_sparsity_ellblock_single_complex(A, threshold);
            break;
        case double_complex:
            return bml_get_sparsity_ellblock_double_complex(A, threshold);
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
