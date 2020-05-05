#include "bml_getters_csr.h"
#include "../bml_introspection.h"
#include "../bml_logger.h"
#include "../bml_allocate.h"
#include "bml_types_csr.h"
#include "../../macros.h"
#include "../../typed.h"

#include <complex.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/** Return a single row element.
 *
 * \param arow The csr row
 * \param j The column index
 * \return The row element
 */
void *TYPED_FUNC(
    csr_get_row_element) (
    csr_sparse_row_t * arow,
    int j)
{
    static REAL_T MINUS_ONE = -1;
    static REAL_T ZERO = 0.;

    int *cols = arow->cols_;
    REAL_T *vals = (REAL_T *) arow->vals_;

    if (j < 0)
    {
        LOG_ERROR("column index is < 0\n");
        return &MINUS_ONE;
    }
    for (int pos = 0; pos < arow->NNZ_; pos++)
    {
        if (cols[pos] == j)
        {
            return &vals[pos];
        }
    }
    return &ZERO;
}

/** Return a single matrix element.
 *
 * \param A The bml matrix
 * \param i The row index
 * \param j The column index
 * \return The matrix element
 */
void *TYPED_FUNC(
    bml_get_csr) (
    bml_matrix_csr_t * A,
    int i,
    int j)
{
    static REAL_T MINUS_ONE = -1;

    if (i < 0 || i >= A->N_)
    {
        LOG_ERROR("row index out of bounds\n");
        return &MINUS_ONE;
    }
    return TYPED_FUNC(csr_get_row_element) (A->data_[i], j);
}

/** Get the column indexes of a matrix row.
 *
 *  \ingroup getters
 *
 *  \param arow The matrix row
 *
 */
int *TYPED_FUNC(
    csr_get_column_indexes) (
    csr_sparse_row_t * arow)
{
    return arow->cols_;
}

/** Get the column entries of a matrix row.
 *
 *  \ingroup getters
 *
 *  \param arow The matrix row
 *
 */
void *TYPED_FUNC(
    csr_get_column_entries) (
    csr_sparse_row_t * arow)
{
    return arow->vals_;
}

/** Get the number of nonzeros of matrix row.
 *
 *  \ingroup getters
 *
 *  \param arow The matrix row
 *
 */
int TYPED_FUNC(
    csr_get_nnz) (
    csr_sparse_row_t * arow)
{
    return arow->NNZ_;
}

/** Get a sparse row i of matrix A.
 *
 *  \ingroup getters
 *
 *  \param A The matrix which takes row i
 *  \param i The index of the row to get
 *  \param cols Array of column indexes
 *  \param vals Array of row values
 *  \param nnz number of nonzero entries
 *
 */
void TYPED_FUNC(
    bml_get_sparse_row_csr) (
    bml_matrix_csr_t * A,
    int i,
    int **cols,
    void **vals,
    int *nnz)
{
    int nz = TYPED_FUNC(csr_get_nnz) (A->data_[i]);
    int *colids = malloc(nz * sizeof(int));
    REAL_T *colvals = malloc(nz * sizeof(REAL_T));
    memcpy(colids, (int *) TYPED_FUNC(csr_get_column_indexes) (A->data_[i]),
           nz * sizeof(int));
    memcpy(colvals,
           (REAL_T *) TYPED_FUNC(csr_get_column_entries) (A->data_[i]),
           nz * sizeof(REAL_T));

    *cols = colids;
    *vals = colvals;
    *nnz = nz;
}

/** Get a dense row i of matrix A. Assumes A is a square matrix.
 *
 *  \ingroup getters
 *
 *  \param A The matrix which takes row i
 *  \param i The index of the row to get
 *  \param row Array to copy row
 *
 */
void *TYPED_FUNC(
    bml_get_row_csr) (
    bml_matrix_csr_t * A,
    int i)
{
    int nnz = TYPED_FUNC(csr_get_nnz) (A->data_[i]);
    int A_N = A->N_;
    int *cols = TYPED_FUNC(csr_get_column_indexes) (A->data_[i]);
    REAL_T *vals = TYPED_FUNC(csr_get_column_entries) (A->data_[i]);

    REAL_T *row = bml_allocate_memory(sizeof(REAL_T) * A_N);

    // loop over column entries to copy row data
    for (int j = 0; j < nnz; j++)
    {
        row[cols[j]] = vals[j];
    }

    return row;
}

/** Get the diagonal of matrix A.
 *
 *  \ingroup getters
 *
 *  \param A The matrix which takes row i
 *  \param Diagonal Array to copy the diagonal
 *
 */
void *TYPED_FUNC(
    bml_get_diagonal_csr) (
    bml_matrix_csr_t * A)
{
    int A_N = A->N_;
    REAL_T *diagonal = calloc(A_N, sizeof(REAL_T));

    for (int i = 0; i < A_N; i++)
    {
        diagonal[i] =
            *((REAL_T *) TYPED_FUNC(csr_get_row_element) (A->data_[i], i));
    }
    return diagonal;
}
