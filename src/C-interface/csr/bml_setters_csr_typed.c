#include "../../macros.h"
#include "../../typed.h"
#include "../bml_introspection.h"
#include "../bml_allocate.h"
#include "../bml_types.h"
#include "bml_setters_csr.h"
#include "bml_types_csr.h"

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

/** Set new row entry i,j (asuming there is no existing entry i,j in the row)
 *
 *  \ingroup setters
 *
 *  \param arow The row which takes entry j
 *  \param j The column index
 *  \param element The element to be added
 *  \WARNING sets an element from scratch
 *  \todo set element new.
 *
 *
 */
void TYPED_FUNC(
    csr_set_row_element_new) (
    csr_sparse_row_t * arow,
    const int j,
    const void *element)
{
    const int pos = arow->NNZ_;
    int *cols = arow->cols_;
    REAL_T *vals = (REAL_T *) arow->vals_;
    // increment nnz counter and reallocate memory if needed
    arow->NNZ_++;
    if (arow->NNZ_ > arow->alloc_size_)
    {
        arow->alloc_size_ *= EXPAND_FACT;
        arow->cols_ =
            bml_reallocate_memory(cols, sizeof(int) * arow->alloc_size_);
        arow->vals_ =
            bml_reallocate_memory(vals, sizeof(REAL_T) * arow->alloc_size_);
        cols = arow->cols_;
        vals = (REAL_T *) arow->vals_;
    }
    cols[pos] = j;
    vals[pos] = *((REAL_T *) element);
}

/** Set row entry i,j. Overwrites existing entry
 *
 *  \ingroup setters
 *
 *  \param arow The row which takes entry j
 *  \param j The column index
 *  \param element The element to be added
 *  \WARNING sets an element from scratch.
 *
 *
 */
void TYPED_FUNC(
    csr_set_row_element) (
    csr_sparse_row_t * arow,
    const int j,
    const void *element)
{
    const int annz = arow->NNZ_;
    int *cols = arow->cols_;
    REAL_T *vals = (REAL_T *) arow->vals_;
    int found = 0;
    // check to see if column entry exists
    for (int l = 0; l < annz; l++)
    {
        if (cols[l] == j)
        {
            vals[l] = *((REAL_T *) element);
            found = 1;
            break;
        }
    }
    if (!found)
    {
        TYPED_FUNC(csr_set_row_element_new) (arow, j, element);
    }
}

/** Set row entries from dense array
 *
 *  \ingroup setters
 *
 *  \param arow The row to be set
 *  \param count The number of entries to set
 *  \param cols The column indexes
 *  \param vals The row entries
 *  \WARNING sets a row from scratch
 *
 *
 */
void TYPED_FUNC(
    csr_set_row) (
    csr_sparse_row_t * arow,
    const int count,
    void *rowvals,
    const double threshold)
{
    int *index = arow->cols_;
    REAL_T *vals = rowvals;
    REAL_T *data = (REAL_T *) arow->vals_;
    // reallocate memory if needed
    if (count > arow->alloc_size_)
    {
        arow->alloc_size_ = count;
        arow->cols_ =
            bml_reallocate_memory(index, sizeof(int) * arow->alloc_size_);
        arow->vals_ =
            bml_reallocate_memory(data, sizeof(REAL_T) * arow->alloc_size_);
        index = arow->cols_;
        data = (REAL_T *) arow->vals_;
    }
    // set entries
    arow->NNZ_ = 0;
    for (int j = 0; j < count; j++)
    {
        if (ABS(vals[j]) > threshold)
        {
            index[arow->NNZ_] = j;
            data[arow->NNZ_++] = vals[j];
        }
    }
}

/** Set (new) element i,j asuming there's no resetting of any element of A.
 *
 *  \ingroup setters
 *
 *  \param A The matrix which takes row i
 *  \param i The column index
 *  \param j The row index
 *  \param value The element to be added
 *  \WARNING sets an element from scratch
 *
 *
 */
void TYPED_FUNC(
    bml_set_element_new_csr) (
    bml_matrix_csr_t * A,
    const int i,
    const int j,
    void *element)
{
    // Insert new entry into row i.
    // Use the pointer to row i directly, since there
    // may be reallocation of memory
    TYPED_FUNC(csr_set_row_element_new) (A->data_[i], j, element);
}


/** Set element i,j of matrix A.
 *
 *  \ingroup setters
 *
 *  \param A The matrix which takes row i
 *  \param i The column index
 *  \param j The row index
 *  \param value The element to be set
 *  \WARNING sets an element from scratch
 *
 *
 */
void TYPED_FUNC(
    bml_set_element_csr) (
    bml_matrix_csr_t * A,
    const int i,
    const int j,
    void *element)
{

    // Insert new entry into row i.
    // Use the pointer to row i directly, since there
    // may be reallocation of memory
    TYPED_FUNC(csr_set_row_element) (A->data_[i], j, element);

}

/** Set row i of matrix A. Assume A is a square matrix.
 *
 *  \ingroup setters
 *
 *  \param A The matrix which takes row i
 *  \param i The index of the row to be set
 *  \param row The row to be set
 *  \param threshold The threshold value to be set
 *
 */
void TYPED_FUNC(
    bml_set_row_csr) (
    bml_matrix_csr_t * A,
    const int i,
    void *rowvals,
    const double threshold)
{
    const int A_N = A->N_;
    csr_sparse_row_t *arow = A->data_[i];
    // reset nnz row count to zero (in case row is not empty)
    arow->NNZ_ = 0;

    // set row values
    TYPED_FUNC(csr_set_row) (arow, A_N, rowvals, threshold);
}

/** Set diagonal of matrix A.
 *
 *  \ingroup setters
 *
 *  \param A The matrix which takes diag
 *  \param diag The diagonal to be set
 *  \param threshold The threshold value to be used
 */
void TYPED_FUNC(
    bml_set_diagonal_csr) (
    bml_matrix_csr_t * A,
    void *diag,
    const double threshold)
{
    REAL_T *diagonal = diag;
    int A_N = A->N_;

    // loop over rows
    for (int i = 0; i < A_N; i++)
    {
        REAL_T dval = ABS(diagonal[i]) > threshold ? diagonal[i] : 0.0;
        TYPED_FUNC(csr_set_row_element) (A->data_[i], i, &dval);
    }
}

/** Set (sparse) row entries
 *
 *  \ingroup setters
 *
 *  \param arow The row to be set
 *  \param count The number of entries to set
 *  \param cols The column indexes
 *  \param vals The row entries
 *  \WARNING sets a row from scratch
 *
 *
 */
void TYPED_FUNC(
    csr_set_sparse_row) (
    csr_sparse_row_t * arow,
    const int count,
    const int *cols,
    void *rowvals,
    const double threshold)
{
    int *index = arow->cols_;
    REAL_T *vals = rowvals;
    REAL_T *data = (REAL_T *) arow->vals_;
    // reallocate memory if needed
    if (count > arow->alloc_size_)
    {
        arow->alloc_size_ = count;
        arow->cols_ =
            bml_reallocate_memory(index, sizeof(int) * arow->alloc_size_);
        arow->vals_ =
            bml_reallocate_memory(data, sizeof(REAL_T) * arow->alloc_size_);
        index = arow->cols_;
        data = (REAL_T *) arow->vals_;
    }
    // set entries
    arow->NNZ_ = 0;
    for (int j = 0; j < count; j++)
    {
        if (ABS(vals[j]) > threshold)
        {
            index[arow->NNZ_] = cols[j];
            data[arow->NNZ_++] = vals[j];
        }
    }
}


/** Set row i of matrix A.
 *
 *  \ingroup setters
 *
 *  \param A The matrix which takes row i
 *  \param i The index of the row to be set
 *  \param count The number of entries to be set
 *  \param cols The column index data to be set
 *  \param vals The row data to be set
 *  \param threshold The threshold value to be set
 *
 */
void TYPED_FUNC(
    bml_set_sparse_row_csr) (
    bml_matrix_csr_t * A,
    const int i,
    const int count,
    const int *cols,
    void *vals,
    const double threshold)
{

    // set row entries
    TYPED_FUNC(csr_set_sparse_row) (A->data_[i], count, cols, vals,
                                    threshold);

}
