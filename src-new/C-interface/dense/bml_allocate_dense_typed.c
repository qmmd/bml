#include "../typed.h"
#include "bml_allocate.h"
#include "bml_allocate_dense.h"
#include "bml_types.h"
#include "bml_types_dense.h"

/** Allocate the zero matrix.
 *
 *  Note that the matrix \f$ a \f$ will be newly allocated. If it is
 *  already allocated then the matrix will be deallocated in the
 *  process.
 *
 *  \ingroup allocate_group
 *
 *  \param N The matrix size.
 *  \return The matrix.
 */
bml_matrix_dense_t *TYPED_FUNC(
    bml_zero_matrix_dense) (
    const int N)
{
    bml_matrix_dense_t *A = NULL;

    A = bml_allocate_memory(sizeof(bml_matrix_dense_t));
    A->matrix_type = dense;
    A->matrix_precision = MATRIX_PRECISION;
    A->N = N;
    A->matrix = bml_allocate_memory(sizeof(REAL_T) * N * N);
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
 *  \return The matrix.
 */
bml_matrix_dense_t *TYPED_FUNC(
    bml_random_matrix_dense) (
    const int N)
{
    bml_matrix_dense_t *A = TYPED_FUNC(bml_zero_matrix_dense) (N);
    REAL_T *A_dense = NULL;

    A_dense = A->matrix;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            A_dense[i + N * j] = rand() / (REAL_T) RAND_MAX;
        }
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
 *  \return The matrix.
 */
bml_matrix_dense_t *TYPED_FUNC(
    bml_identity_matrix_dense) (
    const int N)
{
    bml_matrix_dense_t *A = TYPED_FUNC(bml_zero_matrix_dense) (N);
    REAL_T *A_dense = NULL;

    A_dense = A->matrix;
    for (int i = 0; i < N; i++)
    {
        A_dense[i + N * i] = 1;
    }
    return A;
}