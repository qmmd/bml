#include "../../macros.h"
#include "../typed.h"
#include "bml_allocate.h"
#include "bml_allocate_dense.h"
#include "bml_types.h"
#include "bml_types_dense.h"

#include <complex.h>
#include <string.h>

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
    memset(A->matrix, 0.0, A->N * A->N * sizeof(REAL_T));
}

/** Allocate the zero matrix.
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
    bml_zero_matrix_dense) (
    const int N,
    const bml_distribution_mode_t distrib_mode)
{
    bml_matrix_dense_t *A = NULL;

    A = bml_allocate_memory(sizeof(bml_matrix_dense_t));
    A->matrix_type = dense;
    A->matrix_precision = MATRIX_PRECISION;
    A->N = N;
    A->distribution_mode = distrib_mode;
    A->matrix = bml_allocate_memory(sizeof(REAL_T) * N * N);
    A->domain = bml_default_domain(N, N, distrib_mode);
    A->domain2 = bml_default_domain(N, N, distrib_mode);
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
    const int N,
    const int M,
    const bml_distribution_mode_t distrib_mode)
{
    bml_matrix_dense_t *A =
        TYPED_FUNC(bml_zero_matrix_dense) (N, distrib_mode);
    REAL_T *A_dense = A->matrix;
#pragma omp parallel for default(none) shared(A_dense)
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
 */
bml_matrix_dense_t *TYPED_FUNC(
    bml_random_matrix_dense) (
    const int N,
    const bml_distribution_mode_t distrib_mode)
{
    bml_matrix_dense_t *A =
        TYPED_FUNC(bml_zero_matrix_dense) (N, distrib_mode);
    REAL_T *A_dense = A->matrix;
#pragma omp parallel for default(none) shared(A_dense)
    for (int i = 0; i < N * N; i++)
    {
        A_dense[i] = rand() / (double) RAND_MAX;
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
 *  \param N The matrix size.
 *  \param distrib_mode The distribution mode.
 *  \return The matrix.
 */
bml_matrix_dense_t *TYPED_FUNC(
    bml_identity_matrix_dense) (
    const int N,
    const bml_distribution_mode_t distrib_mode)
{
    bml_matrix_dense_t *A =
        TYPED_FUNC(bml_zero_matrix_dense) (N, distrib_mode);
    REAL_T *A_dense = A->matrix;
#pragma omp parallel for default(none) shared(A_dense)
    for (int i = 0; i < N; i++)
    {
        A_dense[ROWMAJOR(i, i, N, N)] = 1;
    }
    return A;
}
