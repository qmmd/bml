#include "bml_allocate.h"
#include "bml_introspection.h"
#include "bml_logger.h"
#include "bml_parallel.h"
#include "dense/bml_allocate_dense.h"
#include "ellpack/bml_allocate_ellpack.h"
#include "ellblock/bml_allocate_ellblock.h"
#include "ellsort/bml_allocate_ellsort.h"
#include "csr/bml_allocate_csr.h"
#ifdef DO_MPI
#include "distributed2d/bml_allocate_distributed2d.h"
#endif

#include <errno.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/** Check if matrix is allocated.
 *
 * \ingroup allocate_group_C
 *
 * \param A[in,out] Matrix
 * \return \f$ > 0 \f$ if allocated, else -1
 */
int
bml_allocated(
    bml_matrix_t * A)
{
    return bml_get_N(A);
}

/** Allocate and zero a chunk of memory.
 *
 * \ingroup allocate_group_C
 *
 * \param size The size of the memory.
 * \return A pointer to the allocated chunk.
 */
void *
bml_allocate_memory(
    size_t size)
{
#if defined(INTEL_OPT)
    char *ptr = _mm_malloc(size, MALLOC_ALIGNMENT);
#pragma omp parallel for simd
#pragma vector aligned
    for (size_t i = 0; i < size; i++)
    {
        __assume_aligned(ptr, MALLOC_ALIGNMENT);
        ptr[i] = 0;
    }
#elif defined(HAVE_POSIX_MEMALIGN)
    char *ptr;
    posix_memalign((void **) &ptr, MALLOC_ALIGNMENT, size);
    for (size_t i = 0; i < size; i++)
    {
        ptr[i] = 0;
    }
#else
    void *ptr = calloc(1, size);
#endif

    if (ptr == NULL)
    {
        LOG_ERROR("error allocating memory of size %d: %s\n", size,
                  strerror(errno));
    }
    return (void *) ptr;
}

/** Allocate a chunk of memory without initialization.
 *
 * \ingroup allocate_group_C
 *
 * \param size The size of the memory.
 * \return A pointer to the allocated chunk.
 */
void *
bml_noinit_allocate_memory(
    size_t size)
{
#if defined(INTEL_OPT)
    void *ptr = _mm_malloc(size, MALLOC_ALIGNMENT);
#elif defined(HAVE_POSIX_MEMALIGN)
    void *ptr;
    posix_memalign(&ptr, MALLOC_ALIGNMENT, size);
#else
    void *ptr = malloc(size);
#endif
    if (ptr == NULL)
    {
        LOG_ERROR("error allocating memory: %s\n", strerror(errno));
    }
    return ptr;
}

/** Reallocate a chunk of memory.
 *
 * \ingroup allocate_group_C
 *
 * \param size The size of the memory.
 * \return A pointer to the reallocated chunk.
 */
void *
bml_reallocate_memory(
    void *ptr,
    const size_t size)
{
    void *ptr_new = realloc(ptr, size);
    if (ptr_new == NULL)
    {
        LOG_ERROR("error reallocating memory: %s\n", strerror(errno));
    }
    return ptr_new;
}

/** Deallocate a chunk of memory.
 *
 * \ingroup allocate_group_C
 *
 * \param ptr A pointer to the previously allocated chunk.
 */
void
bml_free_memory(
    void *ptr)
{
#ifdef INTEL_OPT
    _mm_free(ptr);
#else
    free(ptr);
#endif
}

/** De-allocate a chunk of memory that was allocated inside a C
 * function. This is used by the Fortran bml_free_C interface. Note
 * the "pointer to pointer" in the API.
 *
 * \ingroup allocate_group_C
 *
 * \param ptr A pointer to the previously allocated chunk.
 */
void
bml_free_ptr(
    void **ptr)
{
    bml_free_memory(*ptr);
}

/** Deallocate a matrix.
 *
 * \ingroup allocate_group_C
 *
 * \param A[in,out] The matrix.
 */
void
bml_deallocate(
    bml_matrix_t ** A)
{
    if (A == NULL)
    {
        LOG_DEBUG("A is NULL\n");
    }
    else if (*A == NULL)
    {
        LOG_DEBUG("*A is NULL\n");
    }
    else
    {
        LOG_DEBUG("deallocating bml matrix\n");
        switch (bml_get_type(*A))
        {
            case dense:
                bml_deallocate_dense(*A);
                break;
            case ellpack:
                bml_deallocate_ellpack(*A);
                break;
            case ellsort:
                bml_deallocate_ellsort(*A);
                break;
            case ellblock:
                bml_deallocate_ellblock(*A);
                break;
            case csr:
                bml_deallocate_csr(*A);
                break;
#ifdef DO_MPI
            case distributed2d:
                bml_deallocate_distributed2d(*A);
                break;
#endif
            default:
                LOG_ERROR("unknown matrix type (%d)\n", bml_get_type(*A));
                break;
        }
        *A = NULL;
    }
}

/** Deallocate a domain.
 *
 * \ingroup allocate_group_C
 *
 * \param D[in,out] The domain.
 */
void
bml_deallocate_domain(
    bml_domain_t * D)
{
    bml_free_memory(D->localRowMin);
    bml_free_memory(D->localRowMax);
    bml_free_memory(D->localRowExtent);
    bml_free_memory(D->localDispl);
    bml_free_memory(D->localElements);
    bml_free_memory(D);
}

/** Clear a matrix.
 *
 * \ingroup allocate_group_C
 *
 * \param A[in,out] The matrix.
 */
void
bml_clear(
    bml_matrix_t * A)
{
    switch (bml_get_type(A))
    {
        case dense:
            bml_clear_dense(A);
            break;
        case ellpack:
            bml_clear_ellpack(A);
            break;
        case ellsort:
            bml_clear_ellsort(A);
            break;
        case ellblock:
            bml_clear_ellblock(A);
            break;
        case csr:
            bml_clear_csr(A);
            break;
#ifdef DO_MPI
        case distributed2d:
            bml_clear_distributed2d(A);
            break;
#endif
        default:
            LOG_ERROR("unknown matrix type (%d)\n", bml_get_type(A));
            break;
    }
}

/** Allocate a matrix without initializing.
 *
 *  Note that the matrix \f$ A \f$ will be newly allocated. The
 *  function does not check whether the matrix is already allocated.
 *
 *  \ingroup allocate_group_C
 *
 *  \param matrix_type The matrix type.
 *  \param matrix_precision The precision of the matrix.
 *  \param matrix_dimension The matrix size.
 *  \param distrib_mode The distribution mode.
 *  \return The matrix.
 */
bml_matrix_t *
bml_noinit_rectangular_matrix(
    bml_matrix_type_t matrix_type,
    bml_matrix_precision_t matrix_precision,
    bml_matrix_dimension_t matrix_dimension,
    bml_distribution_mode_t distrib_mode)
{
    LOG_DEBUG("noinit matrix of size %d (or zero matrix for dense)\n",
              matrix_dimension.N_rows);
#ifdef DO_MPI
    if (distrib_mode == distributed)
        return bml_zero_matrix_distributed2d(matrix_type, matrix_precision,
                                             matrix_dimension.N_rows,
                                             matrix_dimension.N_nz_max);
    else
#endif
        switch (matrix_type)
        {
            case dense:
                return bml_zero_matrix_dense(matrix_precision,
                                             matrix_dimension, distrib_mode);
                break;
            case ellpack:
                return bml_noinit_matrix_ellpack(matrix_precision,
                                                 matrix_dimension,
                                                 distrib_mode);
                break;
            case ellsort:
                return bml_noinit_matrix_ellsort(matrix_precision,
                                                 matrix_dimension,
                                                 distrib_mode);
                break;
            case ellblock:
                return bml_noinit_matrix_ellblock(matrix_precision,
                                                  matrix_dimension,
                                                  distrib_mode);
                break;
            case csr:
                return bml_noinit_matrix_csr(matrix_precision,
                                             matrix_dimension, distrib_mode);
                break;
            default:
                LOG_ERROR("unknown matrix type\n");
                break;
        }
    return NULL;
}

/** Allocate a block matrix
 *
 * \param matrix_type The matrix type.
 * \param matrix_precision The precision of the matrix.
 * \param NB The number of blocks in a row.
 * \param bsizes The sizes of each block
 * \param MB The number of non-zeroes blocks per row.
 * \param distrib_mode The distribution mode.
 * \return The matrix.
 */
bml_matrix_t *
bml_block_matrix(
    bml_matrix_type_t matrix_type,
    bml_matrix_precision_t matrix_precision,
    int NB,
    int MB,
    int M,
    int *bsizes,
    bml_distribution_mode_t distrib_mode)
{
    LOG_DEBUG("block matrix with %d blocks\n", NB);
    switch (matrix_type)
    {
        case ellpack:
            return bml_block_matrix_ellblock(matrix_precision,
                                             NB, MB, M, bsizes, distrib_mode);
            break;
        case ellblock:
            return bml_block_matrix_ellblock(matrix_precision,
                                             NB, MB, M, bsizes, distrib_mode);
            break;
        default:
            LOG_ERROR("unsupported matrix type (type ID %d)\n", matrix_type);
            break;
    }
    return NULL;
}

/** Allocate a matrix without initializing.
 *
 *  Note that the matrix \f$ A \f$ will be newly allocated. The
 *  function does not check whether the matrix is already allocated.
 *
 *  \ingroup allocate_group_C
 *
 *  \param matrix_type The matrix type.
 *  \param matrix_precision The precision of the matrix.
 *  \param N The matrix size.
 *  \param M The number of non-zeroes per row.
 *  \param distrib_mode The distribution mode.
 *  \return The matrix.
 */
bml_matrix_t *
bml_noinit_matrix(
    bml_matrix_type_t matrix_type,
    bml_matrix_precision_t matrix_precision,
    int N,
    int M,
    bml_distribution_mode_t distrib_mode)
{
    bml_matrix_dimension_t matrix_dimension = { N, N, M };
    return bml_noinit_rectangular_matrix(matrix_type, matrix_precision,
                                         matrix_dimension, distrib_mode);
}

/** Allocate the zero matrix.
 *
 *  Note that the matrix \f$ A \f$ will be newly allocated. The
 *  function does not check whether the matrix is already allocated.
 *
 *  \ingroup allocate_group_C
 *
 *  \param matrix_type The matrix type.
 *  \param matrix_precision The precision of the matrix.
 *  \param N The matrix size.
 *  \param M The number of non-zeroes per row.
 *  \param distrib_mode The distribution mode.
 *  \return The matrix.
 */
bml_matrix_t *
bml_zero_matrix(
    bml_matrix_type_t matrix_type,
    bml_matrix_precision_t matrix_precision,
    int N,
    int M,
    bml_distribution_mode_t distrib_mode)
{
    LOG_DEBUG("zero matrix of size %d\n", N);

#ifdef DO_MPI
    if (distrib_mode == distributed)
        return bml_zero_matrix_distributed2d(matrix_type, matrix_precision, N,
                                             M);
    else
#endif
    {
        bml_matrix_dimension_t matrix_dimension = { N, N, M };
        switch (matrix_type)
        {
            case dense:
                return bml_zero_matrix_dense(matrix_precision,
                                             matrix_dimension, distrib_mode);
                break;
            case ellpack:
                return bml_zero_matrix_ellpack(matrix_precision, N, M,
                                               distrib_mode);
                break;
            case ellsort:
                return bml_zero_matrix_ellsort(matrix_precision, N, M,
                                               distrib_mode);
                break;
            case ellblock:
                return bml_zero_matrix_ellblock(matrix_precision, N, M,
                                                distrib_mode);
                break;
            case csr:
                return bml_zero_matrix_csr(matrix_precision, N, M,
                                           distrib_mode);
                break;
            default:
                LOG_ERROR("unknown matrix type\n");
                break;
        }
    }
    return NULL;
}

/** Allocate a random matrix.
 *
 *  Note that the matrix \f$ A \f$ will be newly allocated. The
 *  function does not check whether the matrix is already allocated.
 *
 *  \ingroup allocate_group_C
 *
 *  \param matrix_type The matrix type.
 *  \param matrix_precision The precision of the matrix.
 *  \param N The matrix size.
 *  \param M The number of non-zeroes per row.
 *  \param distrib_mode The distribution mode.
 *  \return The matrix.
 */
bml_matrix_t *
bml_random_matrix(
    bml_matrix_type_t matrix_type,
    bml_matrix_precision_t matrix_precision,
    int N,
    int M,
    bml_distribution_mode_t distrib_mode)
{
#ifdef DO_MPI
    if (distrib_mode == distributed)
        return bml_random_matrix_distributed2d(matrix_type, matrix_precision,
                                               N, M);
    else
#endif
        switch (matrix_type)
        {
            case dense:
                return bml_random_matrix_dense(matrix_precision, N,
                                               distrib_mode);
                break;
            case ellpack:
                return bml_random_matrix_ellpack(matrix_precision, N, M,
                                                 distrib_mode);
                break;
            case ellsort:
                return bml_random_matrix_ellsort(matrix_precision, N, M,
                                                 distrib_mode);
                break;
            case ellblock:
                return bml_random_matrix_ellblock(matrix_precision, N, M,
                                                  distrib_mode);
                break;
            case csr:
                return bml_random_matrix_csr(matrix_precision, N, M,
                                             distrib_mode);
                break;
            default:
                LOG_ERROR("unknown matrix type (type ID %d)\n", matrix_type);
                break;
        }
    return NULL;
}

/** Allocate a banded matrix.
 *
 *  Note that the matrix \f$ A \f$ will be newly allocated. The
 *  function does not check whether the matrix is already allocated.
 *
 *  \ingroup allocate_group_C
 *
 *  \param matrix_type The matrix type.
 *  \param matrix_precision The precision of the matrix.
 *  \param N The matrix size.
 *  \param M The bandwidth of the matrix.
 *  \param distrib_mode The distribution mode.
 *  \return The matrix.
 */
bml_matrix_t *
bml_banded_matrix(
    bml_matrix_type_t matrix_type,
    bml_matrix_precision_t matrix_precision,
    int N,
    int M,
    bml_distribution_mode_t distrib_mode)
{
    LOG_DEBUG("banded matrix of size %d\n", N);
    switch (matrix_type)
    {
        case dense:
            return bml_banded_matrix_dense(matrix_precision, N, M,
                                           distrib_mode);
            break;
        case ellpack:
            return bml_banded_matrix_ellpack(matrix_precision, N, M,
                                             distrib_mode);
            break;
        case ellsort:
            return bml_banded_matrix_ellsort(matrix_precision, N, M,
                                             distrib_mode);
            break;
        case ellblock:
            return bml_banded_matrix_ellblock(matrix_precision, N, M,
                                              distrib_mode);
            break;
        case csr:
            return bml_banded_matrix_csr(matrix_precision, N, M,
                                         distrib_mode);
            break;
        default:
            LOG_ERROR("unknown matrix type (type ID %d)\n", matrix_type);
            break;
    }
    return NULL;
}

/** Allocate the identity matrix.
 *
 *  Note that the matrix \f$ A \f$ will be newly allocated. The
 *  function does not check whether the matrix is already allocated.
 *
 *  \ingroup allocate_group_C
 *
 *  \param matrix_type The matrix type.
 *  \param matrix_precision The precision of the matrix.
 *  \param N The matrix size.
 *  \param M The number of non-zeroes per row.
 *  \param distrib_mode The distribution mode.
 *  \return The matrix.
 */
bml_matrix_t *
bml_identity_matrix(
    bml_matrix_type_t matrix_type,
    bml_matrix_precision_t matrix_precision,
    int N,
    int M,
    bml_distribution_mode_t distrib_mode)
{
    LOG_DEBUG("identity matrix of size %d\n", N);
#ifdef DO_MPI
    if (distrib_mode == distributed)
        return bml_identity_matrix_distributed2d(matrix_type,
                                                 matrix_precision, N, M);
    else
#endif
        switch (matrix_type)
        {
            case dense:
                return bml_identity_matrix_dense(matrix_precision, N,
                                                 distrib_mode);
                break;
            case ellpack:
                return bml_identity_matrix_ellpack(matrix_precision, N, M,
                                                   distrib_mode);
                break;
            case ellsort:
                return bml_identity_matrix_ellsort(matrix_precision, N, M,
                                                   distrib_mode);
                break;
            case ellblock:
                return bml_identity_matrix_ellblock(matrix_precision, N, M,
                                                    distrib_mode);
                break;
            case csr:
                return bml_identity_matrix_csr(matrix_precision, N, M,
                                               distrib_mode);
                break;
            default:
                LOG_ERROR("unknown matrix type (type ID %d)\n", matrix_type);
                break;
        }
    return NULL;
}

/** Allocate a default domain for a bml matrix.
 *
 * \ingroup allocate_group_C
 *
 *  \param N The number of rows
 *  \param M The number of columns
 *  \param distrib_mode The distribution mode
 *  \return The domain
 */
bml_domain_t *
bml_default_domain(
    int N,
    int M,
    bml_distribution_mode_t distrib_mode)
{
    int avgExtent, nleft;
    int nRanks = bml_getNRanks();

    bml_domain_t *domain = bml_allocate_memory(sizeof(bml_domain_t));

    domain->localRowMin = bml_allocate_memory(nRanks * sizeof(int));
    domain->localRowMax = bml_allocate_memory(nRanks * sizeof(int));
    domain->localRowExtent = bml_allocate_memory(nRanks * sizeof(int));
    domain->localDispl = bml_allocate_memory(nRanks * sizeof(int));
    domain->localElements = bml_allocate_memory(nRanks * sizeof(int));

    domain->totalProcs = nRanks;
    domain->totalRows = N;
    domain->totalCols = M;

    domain->globalRowMin = 0;
    domain->globalRowMax = domain->totalRows;
    domain->globalRowExtent = domain->globalRowMax - domain->globalRowMin;

    switch (distrib_mode)
    {
        case sequential:
        {
            // Default - each rank contains entire matrix, even when running distributed
            for (int i = 0; i < nRanks; i++)
            {
                domain->localRowMin[i] = domain->globalRowMin;
                domain->localRowMax[i] = domain->globalRowMax;
                domain->localRowExtent[i] =
                    domain->localRowMax[i] - domain->localRowMin[i];
                domain->localElements[i] =
                    domain->localRowExtent[i] * domain->totalCols;
                domain->localDispl[i] = 0;
            }

        }
            break;

        case distributed:
        {
            // For completely distributed
            avgExtent = N / nRanks;
            domain->maxLocalExtent = ceil((float) N / (float) nRanks);
            domain->minLocalExtent = avgExtent;

            for (int i = 0; i < nRanks; i++)
            {
                domain->localRowExtent[i] = avgExtent;
            }
            nleft = N - nRanks * avgExtent;
            if (nleft > 0)
            {
                for (int i = 0; i < nleft; i++)
                {
                    domain->localRowExtent[i]++;
                }
            }

            /** For first rank */
            domain->localRowMin[0] = domain->globalRowMin;
            domain->localRowMax[0] = domain->localRowExtent[0];

            /** For middle ranks */
            for (int i = 1; i < (nRanks - 1); i++)
            {
                domain->localRowMin[i] = domain->localRowMax[i - 1];
                domain->localRowMax[i] =
                    domain->localRowMin[i] + domain->localRowExtent[i];
            }

            /** For last rank */
            if (nRanks > 1)
            {
                int last = nRanks - 1;
                domain->localRowMin[last] = domain->localRowMax[last - 1];
                domain->localRowMax[last] =
                    domain->localRowMin[last] + domain->localRowExtent[last];
            }

            /** Number of elements and displacement per rank */
            for (int i = 0; i < nRanks; i++)
            {
                domain->localElements[i] =
                    domain->localRowExtent[i] * domain->totalCols;
                domain->localDispl[i] =
                    (i ==
                     0) ? 0 : domain->localDispl[i - 1] +
                    domain->localElements[i - 1];
            }
        }
            break;

        case graph_distributed:
            LOG_ERROR("graph_distibuted not available\n");
            break;

        default:
            LOG_ERROR("unknown distribution method\n");
            break;
    }

    return domain;
}

/** Update a domain for a bml matrix.
 *
 * \ingroup allocate_group_C
 *
 * \param A Matrix with domain
 * \param localPartMin First part on each rank
 * \param localPartMax Last part on each rank
 * \param nnodesInPart Number of nodes in each part
 */
void
bml_update_domain(
    bml_matrix_t * A,
    int *localPartMin,
    int *localPartMax,
    int *nnodesInPart)
{
    switch (bml_get_type(A))
    {
        case dense:
            bml_update_domain_dense(A, localPartMin, localPartMax,
                                    nnodesInPart);
            break;
        case ellpack:
            bml_update_domain_ellpack(A, localPartMin, localPartMax,
                                      nnodesInPart);
            break;
        case ellsort:
            bml_update_domain_ellsort(A, localPartMin, localPartMax,
                                      nnodesInPart);
            break;
        default:
            LOG_ERROR("unknown matrix type (%d)\n", bml_get_type(A));
            break;
    }
}
