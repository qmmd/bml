#include "../../macros.h"
#include "../../typed.h"
#include "../bml_allocate.h"
#include "../bml_types.h"
#include "bml_allocate_csr.h"
#include "bml_types_csr.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/** Clear a csr matrix row.
 *
 * column indexes and non-zero entries are set to zero
 *
 * \ingroup allocate_group
 *
 * \param A The matrix.
 */
void TYPED_FUNC(
    csr_clear_row) (
    csr_sparse_row_t * row)
{
    memset(row->cols_, 0, row->NNZ_ * sizeof(int));
    memset(row->vals_, 0.0, row->NNZ_ * sizeof(REAL_T));
    row->NNZ_ = 0;
}

/** Clear a matrix.
 *
 * total number of non-zeros, column indexes, and values are set to zero.
 *
 * \ingroup allocate_group
 *
 * \param A The matrix.
 */
void TYPED_FUNC(
    bml_clear_csr) (
    bml_matrix_csr_t * A)
{
    const int n = A->N_;
#pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
        TYPED_FUNC(csr_clear_row) ((A->data_)[i]);
    }
    A->TOTNNZ_ = 0;
}

/** Allocate a matrix row with uninitialized values.
 *
 *  Note that the row \f$ a \f$ will be newly allocated. If it is
 *  already allocated then the row will be deallocated in the
 *  process.
 *
 *  \ingroup allocate_group
 *
 *  \param alloc_size The allocation size.
 *  \return The matrix row.
 */
csr_sparse_row_t *TYPED_FUNC(
    csr_noinit_row) (
    const int alloc_size)
{
    csr_sparse_row_t *arow =
        bml_noinit_allocate_memory(sizeof(csr_sparse_row_t));
    const int size =
        INIT_ROW_SPACE >= alloc_size ? INIT_ROW_SPACE : alloc_size;
    arow->cols_ = bml_noinit_allocate_memory(sizeof(int) * size);
    arow->vals_ = bml_noinit_allocate_memory(sizeof(REAL_T) * size);

    arow->alloc_size_ = size;
    arow->NNZ_ = 0;

    return arow;
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
 *  \param N The matrix size.
 *  \param M An estimate of number non-zeroes per row.
 *  \param distrib_mode The distribution mode.
 *  \return The matrix.
 */
bml_matrix_csr_t *TYPED_FUNC(
    bml_noinit_matrix_csr) (
    bml_matrix_dimension_t matrix_dimension,
    bml_distribution_mode_t distrib_mode)
{
    bml_matrix_csr_t *A =
        bml_noinit_allocate_memory(sizeof(bml_matrix_csr_t));
    A->matrix_type = csr;
    A->matrix_precision = MATRIX_PRECISION;
    A->N_ = matrix_dimension.N_rows;
    A->NZMAX_ = matrix_dimension.N_nz_max;
    A->TOTNNZ_ = 0;
    A->distribution_mode = distrib_mode;
    /** allocate csr row data */
    const int N = A->N_;
    A->data_ = bml_noinit_allocate_memory(sizeof(csr_sparse_row_t *) * N);
#pragma omp parallel for
    for (int i = 0; i < N; i++)
    {
        A->data_[i] = TYPED_FUNC(csr_noinit_row) (A->NZMAX_);
    }
    /** allocate hash table **/
    if (distrib_mode == sequential)
    {
        A->table_ = NULL;
    }
    else
    {
        A->table_ = NULL;
//       A->table_ = csr_noinit_table(N);
    }
    /** end allocate hash table **/
/*
    A->domain = bml_default_domain(A->N_, A->NZMAX_, distrib_mode);
    A->domain2 = bml_default_domain(A->N_, A->NZMAX_, distrib_mode);
*/
    return A;
}

/** Allocate a zero matrix row.
 *
 *  Note that the row \f$ a \f$ will be newly allocated. If it is
 *  already allocated then the row will be deallocated in the
 *  process.
 *
 *  \ingroup allocate_group
 *
 *  \param alloc_size The allocation size.
 *  \return The matrix row.
 */
csr_sparse_row_t *TYPED_FUNC(
    csr_zero_row) (
    const int alloc_size)
{
    csr_sparse_row_t *arow =
        bml_noinit_allocate_memory(sizeof(csr_sparse_row_t));
    const int size =
        INIT_ROW_SPACE >= alloc_size ? INIT_ROW_SPACE : alloc_size;
    arow->cols_ = bml_allocate_memory(sizeof(int) * size);
    arow->vals_ = bml_allocate_memory(sizeof(REAL_T) * size);

    arow->alloc_size_ = size;
    arow->NNZ_ = 0;

    return arow;
}

/** Allocate the zero matrix. (Currently assumes sequential case only)
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
 *  \param M An estimate of number non-zeroes per row.
 *  \param distrib_mode The distribution mode.
 *  \return The matrix.
 */
bml_matrix_csr_t *TYPED_FUNC(
    bml_zero_matrix_csr) (
    int N,
    int M,
    bml_distribution_mode_t distrib_mode)
{
    bml_matrix_csr_t *A = bml_allocate_memory(sizeof(bml_matrix_csr_t));
    A->matrix_type = csr;
    A->matrix_precision = MATRIX_PRECISION;
    A->N_ = N;
    A->NZMAX_ = M;
    A->TOTNNZ_ = 0;
    A->distribution_mode = distrib_mode;
    /** allocate csr row data */
    A->data_ = bml_allocate_memory(sizeof(csr_sparse_row_t *) * N);
#pragma omp parallel for
    for (int i = 0; i < N; i++)
    {
        A->data_[i] = TYPED_FUNC(csr_zero_row) (A->NZMAX_);
    }
    /** allocate hash table.
    No need to insert table values since matrix is zero */
    if (distrib_mode == sequential)
    {
        A->table_ = NULL;
    }
    else
    {
        A->table_ = NULL;
//       A->table_ = csr_noinit_table(N);
    }
/*
    A->domain = bml_default_domain(N, M, distrib_mode);
    A->domain2 = bml_default_domain(N, M, distrib_mode);
*/
    return A;
}

/** Allocate a banded random matrix. (Currently assumes sequential case only)
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
 * \param M The bandwidth (the number of non-zero elements per row).
 *  \param distrib_mode The distribution mode.
 *  \return The matrix.
 */
bml_matrix_csr_t *TYPED_FUNC(
    bml_banded_matrix_csr) (
    int N,
    int M,
    bml_distribution_mode_t distrib_mode)
{
    bml_matrix_csr_t *A =
        TYPED_FUNC(bml_zero_matrix_csr) (N, M, distrib_mode);

#pragma omp parallel for
    for (int i = 0; i < N; i++)
    {
        int jind = 0;
        csr_sparse_row_t *row = A->data_[i];
        int *col_indexes = row->cols_;
        REAL_T *row_vals = row->vals_;
        for (int j = (i - M / 2 >= 0 ? i - M / 2 : 0);
             j < (i - M / 2 + M <= N ? i - M / 2 + M : N); j++)
        {
            col_indexes[jind] = j;
            row_vals[jind] = rand() / (REAL_T) RAND_MAX;
            jind++;
        }
        row->NNZ_ = jind;
    }
    /** initialize hash table */
    if (distrib_mode == sequential)
    {
        A->table_ = NULL;
    }
    else
    {
       /** Insert table values here -- not used*/
        A->table_ = NULL;
    }

    return A;
}

/** Allocate a random matrix. (Currently assumes sequential case only)
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
bml_matrix_csr_t *TYPED_FUNC(
    bml_random_matrix_csr) (
    int N,
    int M,
    bml_distribution_mode_t distrib_mode)
{
    bml_matrix_csr_t *A =
        TYPED_FUNC(bml_zero_matrix_csr) (N, M, distrib_mode);

#pragma omp parallel for
    for (int i = 0; i < N; i++)
    {
        int jind = 0;
        csr_sparse_row_t *row = A->data_[i];
        int *col_indexes = row->cols_;
        REAL_T *row_vals = row->vals_;
        for (int j = 0; j < M; j++)
        {
            col_indexes[jind] = j;
            row_vals[jind] = rand() / (REAL_T) RAND_MAX;
            jind++;
        }
        row->NNZ_ = jind;
    }
    /** initialize hash table */
    if (distrib_mode == sequential)
    {
        A->table_ = NULL;
    }
    else
    {
       /** Insert table values here --not used*/
        A->table_ = NULL;
    }

    return A;
}

/** Allocate the identity matrix.(Currently assumes sequential case only)
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
bml_matrix_csr_t *TYPED_FUNC(
    bml_identity_matrix_csr) (
    int N,
    int M,
    bml_distribution_mode_t distrib_mode)
{
    bml_matrix_csr_t *A =
        TYPED_FUNC(bml_zero_matrix_csr) (N, M, distrib_mode);

#pragma omp parallel for
    for (int i = 0; i < N; i++)
    {
        csr_sparse_row_t *row = A->data_[i];
        int *col_indexes = row->cols_;
        REAL_T *row_vals = row->vals_;
        col_indexes[0] = i;
        row_vals[0] = (REAL_T) 1.0;
        row->NNZ_ = 1;
    }
    /** initialize hash table */
    if (distrib_mode == sequential)
    {
        A->table_ = NULL;
    }
    else
    {
       /** Insert table values here --not used*/
        A->table_ = NULL;
    }
    return A;
}
