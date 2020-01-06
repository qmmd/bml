#include "../bml_allocate.h"
#include "../bml_logger.h"
#include "../bml_parallel.h"
#include "../bml_types.h"
#include "bml_allocate_dense.h"
#include "bml_types_dense.h"

/** Deallocate a matrix.
 *
 * \ingroup allocate_group
 *
 * \param A The matrix.
 */
void
bml_deallocate_dense(
    bml_matrix_dense_t * A)
{
#ifdef BML_USE_MAGMA
    magma_queue_destroy(A->queue);
#endif
    bml_deallocate_domain(A->domain);
    bml_deallocate_domain(A->domain2);
#ifdef BML_USE_MAGMA
    magma_int_t ret = magma_free(A->matrix);
    assert(ret == MAGMA_SUCCESS);
#else
    bml_free_memory(A->matrix);
#endif
    bml_free_memory(A);
}

/** Clear a matrix.
 *
 * \ingroup allocate_group
 *
 * \param A The matrix.
 */
void
bml_clear_dense(
    bml_matrix_dense_t * A)
{
    switch (A->matrix_precision)
    {
        case single_real:
            return bml_clear_dense_single_real(A);
            break;
        case double_real:
            return bml_clear_dense_double_real(A);
            break;
        case single_complex:
            return bml_clear_dense_single_complex(A);
            break;
        case double_complex:
            return bml_clear_dense_double_complex(A);
            break;
        default:
            LOG_ERROR("unknown precision (%d)\n", A->matrix_precision);
            break;
    }
}

/** Allocate the zero matrix.
 *
 *  Note that the matrix \f$ a \f$ will be newly allocated. If it is
 *  already allocated then the matrix will be deallocated in the
 *  process.
 *
 *  \ingroup allocate_group
 *
 *  \param matrix_precision The precision of the matrix. The default
 *  is double precision.
 *  \param matrix_dimension The matrix size.
 *  \param distrib_mode The distribution mode.
 *  \return The matrix.
 */
bml_matrix_dense_t *
bml_zero_matrix_dense(
    bml_matrix_precision_t matrix_precision,
    bml_matrix_dimension_t matrix_dimension,
    bml_distribution_mode_t distrib_mode)
{
    switch (matrix_precision)
    {
        case single_real:
            return bml_zero_matrix_dense_single_real(matrix_dimension,
                                                     distrib_mode);
            break;
        case double_real:
            return bml_zero_matrix_dense_double_real(matrix_dimension,
                                                     distrib_mode);
            break;
        case single_complex:
            return bml_zero_matrix_dense_single_complex(matrix_dimension,
                                                        distrib_mode);
            break;
        case double_complex:
            return bml_zero_matrix_dense_double_complex(matrix_dimension,
                                                        distrib_mode);
            break;
        default:
            LOG_ERROR("unknown precision (%d)\n", matrix_precision);
            break;
    }
    return NULL;
}

/** Allocate a banded matrix.
 *
 * Note that the matrix \f$ a \f$ will be newly allocated. If it is
 * already allocated then the matrix will be deallocated in the
 * process.
 *
 * \ingroup allocate_group
 *
 * \param matrix_precision The precision of the matrix. The default
 * is double precision.
 * \param N The matrix size.
 * \param M The bandwidth.
 * \param distrib_mode The distribution mode.
 * \return The matrix.
 */
bml_matrix_dense_t *
bml_banded_matrix_dense(
    bml_matrix_precision_t matrix_precision,
    int N,
    int M,
    bml_distribution_mode_t distrib_mode)
{
    switch (matrix_precision)
    {
        case single_real:
            return bml_banded_matrix_dense_single_real(N, M, distrib_mode);
            break;
        case double_real:
            return bml_banded_matrix_dense_double_real(N, M, distrib_mode);
            break;
        case single_complex:
            return bml_banded_matrix_dense_single_complex(N, M, distrib_mode);
            break;
        case double_complex:
            return bml_banded_matrix_dense_double_complex(N, M, distrib_mode);
            break;
        default:
            LOG_ERROR("unknown precision (%d)\n", matrix_precision);
            break;
    }
    return NULL;
}

/** Allocate a random matrix.
 *
 *  Note that the matrix \f$ a \f$ will be newly allocated. If it is
 *  already allocated then the matrix will be deallocated in the
 *  process.
 *
 *  \ingroup allocate_group
 *
 *  \param matrix_precision The precision of the matrix. The default
 *  is double precision.
 *  \param N The matrix size.
 *  \param distrib_mode The distribution mode.
 *  \return The matrix.
 */
bml_matrix_dense_t *
bml_random_matrix_dense(
    bml_matrix_precision_t matrix_precision,
    int N,
    bml_distribution_mode_t distrib_mode)
{
    switch (matrix_precision)
    {
        case single_real:
            return bml_random_matrix_dense_single_real(N, distrib_mode);
            break;
        case double_real:
            return bml_random_matrix_dense_double_real(N, distrib_mode);
            break;
        case single_complex:
            return bml_random_matrix_dense_single_complex(N, distrib_mode);
            break;
        case double_complex:
            return bml_random_matrix_dense_double_complex(N, distrib_mode);
            break;
        default:
            LOG_ERROR("unknown precision (%d)\n", matrix_precision);
            break;
    }
    return NULL;
}

/** Allocate the identity matrix.
 *
 *  Note that the matrix \f$ a \f$ will be newly allocated. If it is
 *  already allocated then the matrix will be deallocated in the
 *  process.
 *
 *  \ingroup allocate_group
 *
 *  \param matrix_precision The precision of the matrix. The default
 *  is double precision.
 *  \param N The matrix size.
 *  \param distrib_mode The distribution mode
 *  \return The matrix.
 */
bml_matrix_dense_t *
bml_identity_matrix_dense(
    bml_matrix_precision_t matrix_precision,
    int N,
    bml_distribution_mode_t distrib_mode)
{
    switch (matrix_precision)
    {
        case single_real:
            return bml_identity_matrix_dense_single_real(N, distrib_mode);
            break;
        case double_real:
            return bml_identity_matrix_dense_double_real(N, distrib_mode);
            break;
        case single_complex:
            return bml_identity_matrix_dense_single_complex(N, distrib_mode);
            break;
        case double_complex:
            return bml_identity_matrix_dense_double_complex(N, distrib_mode);
            break;
        default:
            LOG_ERROR("unknown precision (%d)\n", matrix_precision);
            break;
    }
    return NULL;
}

/** Update the dense matrix domain.
 *
 * \ingroup allocate_group
 *
 * \param A Matrix with domain
 * \param localPartMin first part on each rank
 * \param localPartMin last part on each rank
 * \param nnodesInPart number of nodes per part
 */
void
bml_update_domain_dense(
    bml_matrix_dense_t * A,
    int *localPartMin,
    int *localPartMax,
    int *nnodesInPart)
{
    bml_domain_t *A_domain = A->domain;

    int myRank = bml_getMyRank();
    int nprocs = bml_getNRanks();

    for (int i = 0; i < nprocs; i++)
    {
        int rtotal = 0;
        for (int j = localPartMin[i]; j <= localPartMax[i]; j++)
        {
            rtotal += nnodesInPart[j - 1];
        }

        if (i == 0)
            A_domain->localRowMin[0] = A_domain->globalRowMin;
        else
            A_domain->localRowMin[i] = A_domain->localRowMax[i - 1];

        A_domain->localRowMax[i] = A_domain->localRowMin[i] + rtotal;
        A_domain->localRowExtent[i] =
            A_domain->localRowMax[i] - A_domain->localRowMin[i];
        A_domain->localElements[i] =
            A_domain->localRowExtent[i] * A_domain->totalCols;

        if (i == 0)
            A_domain->localDispl[0] = 0;
        else
            A_domain->localDispl[i] =
                A_domain->localDispl[i - 1] + A_domain->localElements[i - 1];
    }

    A_domain->minLocalExtent = A_domain->localRowExtent[0];
    A_domain->maxLocalExtent = A_domain->localRowExtent[0];
    for (int i = 1; i < nprocs; i++)
    {
        if (A_domain->localRowExtent[i] < A_domain->minLocalExtent)
            A_domain->minLocalExtent = A_domain->localRowExtent[i];
        if (A_domain->localRowExtent[i] > A_domain->maxLocalExtent)
            A_domain->maxLocalExtent = A_domain->localRowExtent[i];
    }
}
