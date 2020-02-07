#include "../../macros.h"
#include "../../typed.h"
#include "../bml_allocate.h"
#include "../bml_types.h"
#include "bml_allocate_ellblock.h"
#include "bml_setters_ellblock.h"
#include "bml_types_ellblock.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/** Clear a matrix.
 *
 * Numbers of non-zeroes, indeces, and values are set to zero.
 *
 * \ingroup allocate_group
 *
 * \param A The matrix.
 */
void TYPED_FUNC(
    bml_clear_ellblock) (
    bml_matrix_ellblock_t * A)
{
    for (int ib = 0; ib < A->NB; ib++)
        for (int jp = 0; jp < A->nnzb[ib]; jp++)
        {
            int ind = ROWMAJOR(ib, jp, A->NB, A->MB);
            bml_free_memory(A->ptr_value[ind]);
        }
    memset(A->nnzb, 0, A->NB * sizeof(int));
}

bml_matrix_ellblock_t
    * TYPED_FUNC(bml_noinit_matrix_ellblock) (bml_matrix_dimension_t
                                              matrix_dimension,
                                              bml_distribution_mode_t
                                              distrib_mode)
{
    int N = matrix_dimension.N_rows;

    return TYPED_FUNC(bml_block_matrix_ellblock) (bml_get_nb(), bml_get_mb(),
                                                  bml_get_block_sizes(N, 0),
                                                  distrib_mode);
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
bml_matrix_ellblock_t *TYPED_FUNC(
    bml_block_matrix_ellblock) (
    int NB,
    int MB,
    int *bsize,
    bml_distribution_mode_t distrib_mode)
{
    assert(NB > 0);
    assert(MB > 0);
    for (int ib = 0; ib < NB; ib++)
        assert(bsize[ib] < 1e6);

    bml_matrix_ellblock_t *A =
        bml_allocate_memory(sizeof(bml_matrix_ellblock_t));
    A->matrix_type = ellblock;
    A->matrix_precision = MATRIX_PRECISION;
    A->NB = NB;
    A->MB = MB;
    A->bsize = bml_allocate_memory(sizeof(int) * NB);
    memcpy(A->bsize, bsize, NB * sizeof(int));
    A->distribution_mode = distrib_mode;
    A->indexb = bml_allocate_memory(sizeof(int) * NB * MB);
    A->nnzb = bml_allocate_memory(sizeof(int) * NB);
    A->ptr_value = bml_allocate_memory(sizeof(REAL_T *) * NB * MB);

    int N = 0;
    for (int ib = 0; ib < NB; ib++)
        N += bsize[ib];
    A->N = N;
    //printf("bml_block_matrix_ellblock %d %d\n",NB,MB);
    if (bml_get_mb() == 0)
    {
        bml_set_block_sizes(bsize, NB, MB);
    }
    return A;
}

bml_matrix_ellblock_t *TYPED_FUNC(
    bml_zero_matrix_ellblock) (
    int N,
    int M,
    bml_distribution_mode_t distrib_mode)
{
    //use default block values that could have been reset by user
    int *bsize = bml_get_block_sizes(N, M);
    int mb = bml_get_mb();
    int nb = bml_get_nb();

    bml_matrix_ellblock_t *A =
        TYPED_FUNC(bml_block_matrix_ellblock) (nb, mb, bsize, distrib_mode);

    return A;
}

/** Allocate a banded random matrix.
 *
 *  Note that the matrix \f$ a \f$ will be newly allocated. If it is
 *  already allocated then the matrix will be deallocated in the
 *  process.
 *
 *  \ingroup allocate_group
 *
 *  \param N The matrix size.
 *  \param M The number of non-zeroes per row.
 *  \param distrib_mode The distribution mode.
 *  \return The matrix.
 */
bml_matrix_ellblock_t *TYPED_FUNC(
    bml_banded_matrix_ellblock) (
    int N,
    int M,
    bml_distribution_mode_t distrib_mode)
{
    bml_matrix_ellblock_t *A =
        TYPED_FUNC(bml_zero_matrix_ellblock) (N, M, distrib_mode);

    for (int i = 0; i < N; i++)
    {
        for (int j = (i - M / 2 >= 0 ? i - M / 2 : 0);
             j < (i - M / 2 + M <= N ? i - M / 2 + M : N); j++)
        {
            REAL_T value = rand() / (REAL_T) RAND_MAX;
            bml_set_element_new_ellblock(A, i, j, &value);
        }
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
 *
 *  Note: Do not use OpenMP when setting values for a random matrix,
 *  this makes the operation non-repeatable.
 */
bml_matrix_ellblock_t *TYPED_FUNC(
    bml_random_matrix_ellblock) (
    int N,
    int M,
    bml_distribution_mode_t distrib_mode)
{
    //create matrix
    int *bsize = bml_get_block_sizes(N, M);
    int mb = bml_get_mb();
    int nb = bml_get_nb();

    bml_matrix_ellblock_t *A =
        TYPED_FUNC(bml_block_matrix_ellblock) (nb, mb, bsize, distrib_mode);

    //now fill with random values
    int NB = A->NB;
    int MB = A->MB;

    REAL_T **A_ptr_value = (REAL_T **) A->ptr_value;
    int *A_indexb = A->indexb;
    int *A_nnzb = A->nnzb;

    for (int ib = 0; ib < NB; ib++)
    {
        for (int jp = 0; jp < MB; jp++)
        {
            int ind = ROWMAJOR(ib, jp, NB, MB);
            A_indexb[ind] = jp;
            int jb = A_indexb[ind];

            //allocate storage
            int nelements = bsize[ib] * bsize[jb];
            A_ptr_value[ind]
                = bml_noinit_allocate_memory(nelements * sizeof(REAL_T));

            REAL_T *A_value = A_ptr_value[ind];
            assert(A_value != NULL);
            for (int ii = 0; ii < bsize[ib]; ii++)
                for (int jj = 0; jj < bsize[jb]; jj++)
                {
                    A_value[ROWMAJOR(ii, jj, bsize[ib], bsize[jb])] =
                        rand() / (REAL_T) RAND_MAX;
                }
        }
        A_nnzb[ib] = MB;
    }
    return A;
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
bml_matrix_ellblock_t *TYPED_FUNC(
    bml_identity_matrix_ellblock) (
    int N,
    int M,
    bml_distribution_mode_t distrib_mode)
{
    int *bsize = bml_get_block_sizes(N, M);
    int mb = bml_get_mb();
    int nb = bml_get_nb();
    assert(nb > 0);
    assert(mb > 0);

    bml_matrix_ellblock_t *A =
        TYPED_FUNC(bml_block_matrix_ellblock) (nb, mb, bsize, distrib_mode);

    REAL_T **A_ptr_value = (REAL_T **) A->ptr_value;

    int NB = A->NB;
    int MB = A->MB;

    for (int ib = 0; ib < NB; ib++)
    {
        int nelements = bsize[ib] * bsize[ib];
        int ind = ROWMAJOR(ib, 0, NB, MB);
        A_ptr_value[ind] = bml_allocate_memory(nelements * sizeof(REAL_T));
        REAL_T *A_value = A_ptr_value[ind];
        for (int ii = 0; ii < bsize[ib]; ii++)
        {
            A_value[ROWMAJOR(ii, ii, bsize[ib], bsize[ib])] = 1.0;
        }
        A->nnzb[ib] = 1;
        A->indexb[ind] = ib;
    }

    return A;
}
