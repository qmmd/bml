#include "../../macros.h"
#include "../bml_allocate.h"
#include "../bml_logger.h"
#include "../bml_parallel.h"
#include "../bml_types.h"
#include "bml_allocate_ellblock.h"
#include "bml_types_ellblock.h"

#include <assert.h>
#include <stdio.h>

/*
 * variables visible only in that file
 */
static int *s_default_bsize = NULL;
static int s_default_block_dim = 4;
static int s_nb = 0;
static int s_mb = 0;

/*
 *  * function to set values of static variables visible in that file only
 *   * Should be called only once.
 *    */
void
bml_set_block_sizes(
    int *bsize,
    int nb,
    int mb)
{
    printf("bml_set_block_sizes %d %d\n", nb, mb);
    assert(nb > 0);
    assert(mb > 0);

    /*
     * make sure this function is called only once
     */
    assert(s_nb == 0);

    s_nb = nb;
    s_mb = mb;
    s_default_bsize = malloc(nb * sizeof(int));
    for (int ib = 0; ib < nb; ib++)
    {
        s_default_bsize[ib] = bsize[ib];
    }
}

int *
bml_get_block_sizes(
    int N,
    int M)
{
    if (s_default_bsize == NULL)
    {
        assert(M > 0);

        s_nb = N / s_default_block_dim;
        if (s_nb * s_default_block_dim < N)
            s_nb++;
        printf("s_nb = %d\n", s_nb);

        s_mb = M / s_default_block_dim;
        if (s_mb * s_default_block_dim < M)
            s_mb++;
        printf("s_mb = %d\n", s_mb);

        s_default_bsize = malloc(s_nb * sizeof(int));
        for (int ib = 0; ib < s_nb - 1; ib++)
        {
            s_default_bsize[ib] = s_default_block_dim;
        }
        s_default_bsize[s_nb - 1] = N % s_default_block_dim;
    }
    return s_default_bsize;
}

int
bml_get_mb(
    )
{
    return s_mb;
}

int
bml_get_nb(
    )
{
    return s_nb;
}

int
count_nelements(
    bml_matrix_ellblock_t * A)
{
    int nelements = 0;
    for (int ib = 0; ib < A->NB; ib++)
        for (int jp = 0; jp < A->MB; jp++)
        {
            int ind = ROWMAJOR(ib, jp, A->NB, A->MB);
            int jb = A->indexb[ind];
            nelements = nelements + A->bsize[ib] * A->bsize[jb];
        }
    //printf("nelements=%d\n", nelements);
    return nelements;
}

/** Deallocate a matrix.
 *
 * \ingroup allocate_group
 *
 * \param A The matrix.
 */
void
bml_deallocate_ellblock(
    bml_matrix_ellblock_t * A)
{
    for (int ib = 0; ib < A->NB; ib++)
        for (int jp = 0; jp < A->MB; jp++)
        {
            int ind = ROWMAJOR(ib, jp, A->NB, A->MB);
            bml_free_memory(A->ptr_value[ind]);
        }
    bml_free_memory(A->ptr_value);
    bml_free_memory(A->indexb);
    bml_free_memory(A->nnzb);
    bml_free_memory(A->bsize);
    bml_free_memory(A);
}

/** Clear a matrix.
 *
 * \ingroup allocate_group
 *
 * \param A The matrix.
 */
void
bml_clear_ellblock(
    bml_matrix_ellblock_t * A)
{
    switch (A->matrix_precision)
    {
        case single_real:
            bml_clear_ellblock_single_real(A);
            break;
        case double_real:
            bml_clear_ellblock_double_real(A);
            break;
        case single_complex:
            bml_clear_ellblock_single_complex(A);
            break;
        case double_complex:
            bml_clear_ellblock_double_complex(A);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
}

bml_matrix_ellblock_t *
bml_noinit_matrix_ellblock(
    bml_matrix_precision_t matrix_precision,
    bml_matrix_dimension_t matrix_dimension,
    bml_distribution_mode_t distrib_mode)
{
    bml_matrix_ellblock_t *A = NULL;

    switch (matrix_precision)
    {
        case single_real:
            A = bml_noinit_matrix_ellblock_single_real(matrix_dimension,
                                                       distrib_mode);
            break;
        case double_real:
            A = bml_noinit_matrix_ellblock_double_real(matrix_dimension,
                                                       distrib_mode);
            break;
        case single_complex:
            A = bml_noinit_matrix_ellblock_single_complex(matrix_dimension,
                                                          distrib_mode);
            break;
        case double_complex:
            A = bml_noinit_matrix_ellblock_double_complex(matrix_dimension,
                                                          distrib_mode);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
    return A;
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
 *  \param N The matrix size.
 *  \param M The number of non-zeroes per row.
 *  \param distrib_mode The distribution mode.
 *  \return The matrix.
 */
bml_matrix_ellblock_t *
bml_zero_matrix_ellblock(
    bml_matrix_precision_t matrix_precision,
    int N,
    int M,
    bml_distribution_mode_t distrib_mode)
{
    bml_matrix_ellblock_t *A = NULL;

    switch (matrix_precision)
    {
        case single_real:
            A = bml_zero_matrix_ellblock_single_real(N, M, distrib_mode);
            break;
        case double_real:
            A = bml_zero_matrix_ellblock_double_real(N, M, distrib_mode);
            break;
        case single_complex:
            A = bml_zero_matrix_ellblock_single_complex(N, M, distrib_mode);
            break;
        case double_complex:
            A = bml_zero_matrix_ellblock_double_complex(N, M, distrib_mode);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
    return A;
}

bml_matrix_ellblock_t *
bml_block_matrix_ellblock(
    bml_matrix_precision_t matrix_precision,
    int NB,
    int MB,
    int *bsizes,
    bml_distribution_mode_t distrib_mode)
{
    bml_matrix_ellblock_t *A = NULL;

    switch (matrix_precision)
    {
        case single_real:
            A = bml_block_matrix_ellblock_single_real(NB, MB, bsizes,
                                                      distrib_mode);
            break;
        case double_real:
            A = bml_block_matrix_ellblock_double_real(NB, MB, bsizes,
                                                      distrib_mode);
            break;
        case single_complex:
            A = bml_block_matrix_ellblock_single_complex(NB, MB, bsizes,
                                                         distrib_mode);
            break;
        case double_complex:
            A = bml_block_matrix_ellblock_double_complex(NB, MB, bsizes,
                                                         distrib_mode);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
    return A;
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
 *  \param M The number of non-zeroes per row.
 *  \param distrib_mode The distribution mode.
 *  \return The matrix.
 */
bml_matrix_ellblock_t *
bml_random_matrix_ellblock(
    bml_matrix_precision_t matrix_precision,
    int N,
    int M,
    bml_distribution_mode_t distrib_mode)
{
    switch (matrix_precision)
    {
        case single_real:
            return bml_random_matrix_ellblock_single_real(N, M, distrib_mode);
            break;
        case double_real:
            return bml_random_matrix_ellblock_double_real(N, M, distrib_mode);
            break;
        case single_complex:
            return bml_random_matrix_ellblock_single_complex(N, M,
                                                             distrib_mode);
            break;
        case double_complex:
            return bml_random_matrix_ellblock_double_complex(N, M,
                                                             distrib_mode);
            break;
        default:
            LOG_ERROR("unknown precision\n");
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
 *  \param M The number of non-zeroes per row.
 *  \param distrib_mode The distribution mode.
 *  \return The matrix.
 */
bml_matrix_ellblock_t *
bml_identity_matrix_ellblock(
    bml_matrix_precision_t matrix_precision,
    int N,
    int M,
    bml_distribution_mode_t distrib_mode)
{
    switch (matrix_precision)
    {
        case single_real:
            return bml_identity_matrix_ellblock_single_real(N, M,
                                                            distrib_mode);
            break;
        case double_real:
            return bml_identity_matrix_ellblock_double_real(N, M,
                                                            distrib_mode);
            break;
        case single_complex:
            return bml_identity_matrix_ellblock_single_complex(N, M,
                                                               distrib_mode);
            break;
        case double_complex:
            return bml_identity_matrix_ellblock_double_complex(N, M,
                                                               distrib_mode);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
    return NULL;
}

/** Update the ellblock matrix domain.
 *
 * \ingroup allocate_group
 *
 * \param A Matrix with domain
 * \param localPartMin first part on each rank
 * \param localPartMin last part on each rank
 * \param nnodesInPart number of nodes per part
 */
void
bml_update_domain_ellblock(
    bml_matrix_ellblock_t * A,
    int *localPartMin,
    int *localPartMax,
    int *nnodesInPart)
{
    LOG_ERROR("bml_update_domain_ellblock not implemented\n");
}
