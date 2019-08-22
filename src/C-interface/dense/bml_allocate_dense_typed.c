#include "../../macros.h"
#include "../../typed.h"
#include "../bml_allocate.h"
#include "../bml_logger.h"
#include "../bml_types.h"
#include "bml_allocate_dense.h"
#include "bml_types_dense.h"
#include "bml_utilities_dense.h"

#ifdef BML_USE_MAGMA
#include "magma_v2.h"
#endif

#include <complex.h>
#include <string.h>
#include <stdio.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/** Clear the matrix.
 *
 * All values are zeroed.
 *
 * \ingroup allocate_group
 *
 *  \param A The matrix.
 */
void TYPED_FUNC(
    bml_clear_dense) (
    bml_matrix_dense_t * A)
{
#ifdef BML_USE_MAGMA
    MAGMA_T zero = MAGMACOMPLEX(MAKE) (0., 0.);
    MAGMABLAS(laset) (MagmaFull, A->N, A->N, zero, zero, A->matrix, A->ld,
                      A->queue);
#else
    memset(A->matrix, 0.0, A->N * A->ld * sizeof(REAL_T));
#endif
}

/** Allocate the zero matrix.
 *
 *  Note that the matrix \f$ a \f$ will be newly allocated. If it is
 *  already allocated then the matrix will be deallocated in the
 *  process.
 *
 *  \ingroup allocate_group
 *
 *  \param matrix_dimension The matrix size.
 *  \param distrib_mode The distribution mode.
 *  \return The matrix.
 */
bml_matrix_dense_t *TYPED_FUNC(
    bml_zero_matrix_dense) (
    bml_matrix_dimension_t matrix_dimension,
    bml_distribution_mode_t distrib_mode)
{
    bml_matrix_dense_t *A = NULL;

    A = bml_allocate_memory(sizeof(bml_matrix_dense_t));
    A->matrix_type = dense;
    A->matrix_precision = MATRIX_PRECISION;
    A->N = matrix_dimension.N_rows;
    A->distribution_mode = distrib_mode;
#ifdef BML_USE_MAGMA
    A->ld = magma_roundup(matrix_dimension.N_rows, 32);
    int device;
    magma_getdevice(&device);
    magma_queue_create(device, &A->queue);

    magma_int_t ret = MAGMA(malloc) ((MAGMA_T **) & A->matrix,
                                     A->ld * matrix_dimension.N_rows);
    assert(ret == MAGMA_SUCCESS);

    bml_clear_dense(A);
#else
    A->ld = matrix_dimension.N_rows;
    A->matrix =
        bml_allocate_memory(sizeof(REAL_T) * matrix_dimension.N_rows *
                            matrix_dimension.N_rows);
#endif
    A->domain =
        bml_default_domain(matrix_dimension.N_rows, matrix_dimension.N_rows,
                           distrib_mode);
    A->domain2 =
        bml_default_domain(matrix_dimension.N_rows, matrix_dimension.N_rows,
                           distrib_mode);
    return A;
}

/** Allocate a banded matrix.
 *
 * Note that the matrix \f$ a \f$ will be newly allocated. If it is
 * already allocated then the matrix will be deallocated in the
 * process.
 *
 * \ingroup allocate_group
 *
 * \param N The matrix size.
 * \param M The bandwidth (the number of non-zero elements per row).
 * \param distrib_mode The distribution mode.
 * \return The matrix.
 */
bml_matrix_dense_t *TYPED_FUNC(
    bml_banded_matrix_dense) (
    int N,
    int M,
    bml_distribution_mode_t distrib_mode)
{
    bml_matrix_dimension_t matrix_dimension = { N, N, M };
    bml_matrix_dense_t *A =
        TYPED_FUNC(bml_zero_matrix_dense) (matrix_dimension, distrib_mode);
    REAL_T *A_dense = A->matrix;
#pragma omp parallel for shared(A_dense)
    for (int i = 0; i < N; i++)
    {
        for (int j = (i - M / 2 >= 0 ? i - M / 2 : 0);
             j < (i - M / 2 + M <= N ? i - M / 2 + M : N); j++)
        {
            A_dense[ROWMAJOR(i, j, N, N)] = rand() / (double) RAND_MAX;
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
 *  \param N The matrix size.
 *  \param distrib_mode The distribution mode.
 *  \return The matrix.
 *
 *  Note: Do not use OpenMP when setting values for a random matrix,
 *  this makes the operation non-repeatable.
 */
bml_matrix_dense_t *TYPED_FUNC(
    bml_random_matrix_dense) (
    int N,
    bml_distribution_mode_t distrib_mode)
{
    bml_matrix_dimension_t matrix_dimension = { N, N, N };
    bml_matrix_dense_t *A =
        TYPED_FUNC(bml_zero_matrix_dense) (matrix_dimension, distrib_mode);
#ifdef BML_USE_MAGMA
    MAGMA_T *A_dense = malloc(N * N * sizeof(REAL_T));
#else
    REAL_T *A_dense = A->matrix;
#endif

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
#ifdef BML_USE_MAGMA
            A_dense[ROWMAJOR(i, j, N, N)] =
                MAGMACOMPLEX(MAKE) (rand() / (double) RAND_MAX, 0.);
#else
            A_dense[ROWMAJOR(i, j, N, N)] = rand() / (double) RAND_MAX;
#endif
        }
    }

#ifdef BML_USE_MAGMA
    MAGMA(setmatrix) (N, N, A_dense, N, A->matrix, A->ld, A->queue);
    free(A_dense);
#endif
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
 *  \param N The matrix size.
 *  \param distrib_mode The distribution mode.
 *  \return The matrix.
 */
bml_matrix_dense_t *TYPED_FUNC(
    bml_identity_matrix_dense) (
    int N,
    bml_distribution_mode_t distrib_mode)
{
    bml_matrix_dimension_t matrix_dimension = { N, N, N };
    bml_matrix_dense_t *A =
        TYPED_FUNC(bml_zero_matrix_dense) (matrix_dimension, distrib_mode);
#ifdef BML_USE_MAGMA
    MAGMA_T *A_dense = calloc(N * N, sizeof(REAL_T));
#else
    REAL_T *A_dense = A->matrix;
#endif

#pragma omp parallel for shared(A_dense)
    for (int i = 0; i < N; i++)
    {
#ifdef BML_USE_MAGMA
        A_dense[ROWMAJOR(i, i, N, N)] = MAGMACOMPLEX(MAKE) (1, 0);
#else
        A_dense[ROWMAJOR(i, i, N, N)] = 1;
#endif
    }

#ifdef BML_USE_MAGMA
    MAGMA(setmatrix) (N, N, A_dense, N, A->matrix, A->ld, A->queue);
    free(A_dense);
#endif
    return A;
}
