#include "../../macros.h"
#include "../../typed.h"
#include "../bml_allocate.h"
#include "../bml_types.h"
#include "bml_allocate_ellpack.h"
#include "bml_types_ellpack.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/** Deallocate a matrix.
 *
 * \ingroup allocate_group
 *
 * \param A The matrix.
 */
void TYPED_FUNC(
    bml_deallocate_ellpack) (
    bml_matrix_ellpack_t * A)
{
#ifdef USE_OMP_OFFLOAD
    int N = A->N;
    int M = A->M;

    int *A_nnz = A->nnz;
    int *A_index = A->index;
    REAL_T *A_value = A->value;

#pragma omp target exit data map(delete: A_nnz[:N], A_index[:N*M], A_value[:N*M])
#endif

    bml_deallocate_domain(A->domain);
    bml_deallocate_domain(A->domain2);
    bml_free_memory(A->value);
    bml_free_memory(A->index);
    bml_free_memory(A->nnz);
    bml_free_memory(A);
}

/** Clear a matrix.
 *
 * Numbers of non-zeroes, indeces, and values are set to zero.
 *
 * \ingroup allocate_group
 *
 * \param A The matrix.
 */
void TYPED_FUNC(
    bml_clear_ellpack) (
    bml_matrix_ellpack_t * A)
{
    REAL_T *A_value = A->value;
#if defined (USE_OMP_OFFLOAD)
    int *A_index = A->index;
    int *A_nnz = A->nnz;
    int N = A->N;
    int M = A->M;

#pragma omp target teams distribute parallel for
    for (int i = 0; i < N; i++)
    {
        A_nnz[i] = 0;
    }

#pragma omp target teams distribute parallel for collapse(2) schedule (static, 1)
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
            A_index[ROWMAJOR(i, j, N, M)] = 0;
            A_value[ROWMAJOR(i, j, N, M)] = 0.0;
        }
    }
#else // conditional for offload

#ifdef INTEL_OPT
#pragma omp parallel for simd
#pragma vector aligned
    for (int i = 0; i < (A->N * A->M); i++)
    {
        __assume_aligned(A->index, MALLOC_ALIGNMENT);
        __assume_aligned(A_value, MALLOC_ALIGNMENT);
        A->index[i] = 0;
        A_value[i] = 0.0;
    }

#pragma omp parallel for simd
#pragma vector aligned
    for (int i = 0; i < A->N; i++)
    {
        __assume_aligned(A->nnz, MALLOC_ALIGNMENT);
        A->nnz[i] = 0;
    }
#else
    memset(A->nnz, 0, A->N * sizeof(int));
    memset(A->index, 0, A->N * A->M * sizeof(int));
    memset(A->value, 0.0, A->N * A->M * sizeof(REAL_T));
#endif

#endif // conditional for offload
}

/** Allocate a matrix with uninitialized values.
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
bml_matrix_ellpack_t
    * TYPED_FUNC(bml_noinit_matrix_ellpack) (bml_matrix_dimension_t
                                             matrix_dimension,
                                             bml_distribution_mode_t
                                             distrib_mode)
{
    bml_matrix_ellpack_t *A =
        bml_noinit_allocate_memory(sizeof(bml_matrix_ellpack_t));
    A->matrix_type = ellpack;
    A->matrix_precision = MATRIX_PRECISION;
    A->N = matrix_dimension.N_rows;
    A->M = matrix_dimension.N_nz_max;
    A->distribution_mode = distrib_mode;
    A->index = bml_noinit_allocate_memory(sizeof(int) * A->N * A->M);
    A->nnz = bml_allocate_memory(sizeof(int) * A->N);
    A->value = bml_noinit_allocate_memory(sizeof(REAL_T) * A->N * A->M);
    A->domain = bml_default_domain(A->N, A->M, distrib_mode);
    A->domain2 = bml_default_domain(A->N, A->M, distrib_mode);

#if defined(USE_OMP_OFFLOAD)
    int N = A->N;
    int M = A->M;
    int *A_index = A->index;
    int *A_nnz = A->nnz;
    REAL_T *A_value = A->value;

#pragma omp target enter data map(alloc:A_value[:N*M], A_index[:N*M], A_nnz[:N])
#pragma omp target update to(A_value[:N*M], A_index[:N*M], A_nnz[:N])
#endif

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
bml_matrix_ellpack_t *TYPED_FUNC(
    bml_zero_matrix_ellpack) (
    int N,
    int M,
    bml_distribution_mode_t distrib_mode)
{
    assert(M > 0);

    bml_matrix_ellpack_t *A =
        bml_allocate_memory(sizeof(bml_matrix_ellpack_t));
    A->matrix_type = ellpack;
    A->matrix_precision = MATRIX_PRECISION;
    A->N = N;
    A->M = M;
    A->distribution_mode = distrib_mode;
    // need to keep these allocates for host copy
    A->index = bml_allocate_memory(sizeof(int) * N * M);
    A->nnz = bml_allocate_memory(sizeof(int) * N);
    A->value = bml_allocate_memory(sizeof(REAL_T) * N * M);

    REAL_T *A_value = A->value;

    A->domain = bml_default_domain(N, M, distrib_mode);
    A->domain2 = bml_default_domain(N, M, distrib_mode);

#if defined(USE_OMP_OFFLOAD)
    int *A_nnz = A->nnz;
    int *A_index = A->index;
    int NM = N * M;

#pragma omp target enter data map(alloc:A_value[:N*M], A_index[:N*M], A_nnz[:N])

#pragma omp target teams distribute parallel for schedule (static, 1)
    for (int i = 0; i < N; i++)
    {
        A_nnz[i] = 0;
    }

#pragma omp target teams distribute parallel for collapse(2) schedule (static, 1)
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
            A_index[ROWMAJOR(i, j, N, M)] = 0;
            A_value[ROWMAJOR(i, j, N, M)] = 0.0;
        }
    }
#endif
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
 *  \param matrix_precision The precision of the matrix. The default
 *  is double precision.
 *  \param N The matrix size.
 *  \param M The number of non-zeroes per row.
 *  \param distrib_mode The distribution mode.
 *  \return The matrix.
 */
bml_matrix_ellpack_t *TYPED_FUNC(
    bml_banded_matrix_ellpack) (
    int N,
    int M,
    bml_distribution_mode_t distrib_mode)
{
    bml_matrix_ellpack_t *A =
        TYPED_FUNC(bml_zero_matrix_ellpack) (N, M, distrib_mode);

    REAL_T *A_value = A->value;
    int *A_index = A->index;
    int *A_nnz = A->nnz;
    const REAL_T INV_RAND_MAX = 1.0 / (REAL_T) RAND_MAX;

#pragma omp parallel for shared(A_value, A_index, A_nnz)
    for (int i = 0; i < N; i++)
    {
        int jind = 0;
        for (int j = (i - M / 2 >= 0 ? i - M / 2 : 0);
             j < (i - M / 2 + M <= N ? i - M / 2 + M : N); j++)
        {
            A_value[ROWMAJOR(i, jind, N, M)] = rand() * INV_RAND_MAX;
            A_index[ROWMAJOR(i, jind, N, M)] = j;
            jind++;
        }
        A_nnz[i] = jind;
    }
#if defined(USE_OMP_OFFLOAD)
#pragma omp target update to(A_value[:N*M], A_index[:N*M], A_nnz[:N])
#endif
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
bml_matrix_ellpack_t *TYPED_FUNC(
    bml_random_matrix_ellpack) (
    int N,
    int M,
    bml_distribution_mode_t distrib_mode)
{
    bml_matrix_ellpack_t *A =
        TYPED_FUNC(bml_zero_matrix_ellpack) (N, M, distrib_mode);

    REAL_T *A_value = A->value;
    int *A_index = A->index;
    int *A_nnz = A->nnz;
    const REAL_T INV_RAND_MAX = 1.0 / (REAL_T) RAND_MAX;
    for (int i = 0; i < N; i++)
    {
        int jind = 0;
        for (int j = 0; j < M; j++)
        {
            A_value[ROWMAJOR(i, jind, N, M)] = rand() * INV_RAND_MAX;
            A_index[ROWMAJOR(i, jind, N, M)] = j;
            jind++;
        }
        A_nnz[i] = jind;
    }
#if defined(USE_OMP_OFFLOAD)
#pragma omp target update to(A_value[:N*M], A_index[:N*M], A_nnz[:N])
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
 *  \param matrix_precision The precision of the matrix. The default
 *  is double precision.
 *  \param N The matrix size.
 *  \param M The number of non-zeroes per row.
 *  \param distrib_mode The distribution mode.
 *  \return The matrix.
 */
bml_matrix_ellpack_t *TYPED_FUNC(
    bml_identity_matrix_ellpack) (
    int N,
    int M,
    bml_distribution_mode_t distrib_mode)
{
    bml_matrix_ellpack_t *A =
        TYPED_FUNC(bml_zero_matrix_ellpack) (N, M, distrib_mode);

    REAL_T *A_value = A->value;
    int *A_index = A->index;
    int *A_nnz = A->nnz;

#pragma omp parallel for shared(A_value, A_index, A_nnz)
    for (int i = 0; i < N; i++)
    {
#ifdef INTEL_OPT
        __assume_aligned(A_value, MALLOC_ALIGNMENT);
        __assume_aligned(A_index, MALLOC_ALIGNMENT);
        __assume_aligned(A_nnz, MALLOC_ALIGNMENT);
#endif
        A_value[ROWMAJOR(i, 0, N, M)] = (REAL_T) 1.0;
        A_index[ROWMAJOR(i, 0, N, M)] = i;
        A_nnz[i] = 1;
    }
#if defined(USE_OMP_OFFLOAD)
#pragma omp target update to(A_value[:N*M], A_index[:N*M], A_nnz[:N])
#endif
    return A;
}
