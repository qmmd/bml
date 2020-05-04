#include "../../macros.h"
#include "../../typed.h"
#include "../bml_allocate.h"
#include "../bml_parallel.h"
#include "../bml_transpose.h"
#include "../bml_types.h"
#include "bml_allocate_csr.h"
#include "bml_transpose_csr.h"
#include "bml_types_csr.h"
#include "bml_setters_csr.h"
#include "../bml_logger.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/** Transpose a matrix.
 *
 *  \ingroup transpose_group
 *
 *  \param A The matrix to be transposed
 *  \return the transposed A
 */
bml_matrix_csr_t *TYPED_FUNC(
    bml_transpose_new_csr) (
    bml_matrix_csr_t * A)
{
    bml_matrix_dimension_t matrix_dimension = { A->N_, A->N_, A->NZMAX_ };

    bml_matrix_csr_t *B = TYPED_FUNC(bml_noinit_matrix_csr)
        (matrix_dimension, A->distribution_mode);

#ifdef _OPENMP
    omp_lock_t *row_lock =
        (omp_lock_t *) malloc(sizeof(omp_lock_t) * matrix_dimension.N_rows);

#pragma omp parallel for
    for (int i = 0; i < matrix_dimension.N_rows; i++)
    {
        omp_init_lock(&row_lock[i]);
    }
#endif

#pragma omp parallel for                                                \
  shared(matrix_dimension, row_lock)
    for (int i = 0; i < matrix_dimension.N_rows; i++)
    {
        int *cols = A->data_[i]->cols_;
        REAL_T *vals = (REAL_T *) A->data_[i]->vals_;
        const int annz = A->data_[i]->NNZ_;
        for (int pos = 0; pos < annz; pos++)
        {
            const int j = cols[pos];
#ifdef _OPENMP
            omp_set_lock(&row_lock[j]);
#endif
            bml_set_element_new_csr(B, j, i, &vals[pos]);
#ifdef _OPENMP
            omp_unset_lock(&row_lock[j]);
#endif
        }
    }
#ifdef _OPENMP
#pragma omp parallel for
    for (int i = 0; i < matrix_dimension.N_rows; i++)
    {
        omp_destroy_lock(&row_lock[i]);
    }

    free(row_lock);
#endif

    return B;
}

/** swap row entries in position ipos and jpos.
 *
 * column indexes and non-zero entries are swapped
 *
 * \ingroup transpose_group
 *
 * \param A The matrix.
 */
void TYPED_FUNC(
    csr_swap_row_entries) (
    csr_sparse_row_t * row,
    const int ipos,
    const int jpos)
{
    if (ipos == jpos)
        return;

    REAL_T *vals = (REAL_T *) row->vals_;
    int *cols = row->cols_;
    REAL_T tmp = vals[ipos];
    int itmp = cols[ipos];
    /* swap */
    vals[ipos] = vals[jpos];
    vals[jpos] = tmp;
    cols[ipos] = cols[jpos];
    cols[jpos] = itmp;
}

/** Transpose a matrix in place.
 *
 *  \ingroup transpose_group
 *
 *  \param A The matrix to be transposeed
 *  \return the transposed A
 */
void TYPED_FUNC(
    bml_transpose_csr) (
    bml_matrix_csr_t * A)
{
    int N = A->N_;
    int nz_t[N];
    memset(nz_t, 0, sizeof(int) * N);

//#pragma omp parallel for shared(N, M, A_value, A_index, A_nnz)
    // symmetric and nonsymmetric contributions from upper triangular part
    for (int i = 0; i < N; i++)
    {
        int *icols = A->data_[i]->cols_;
        REAL_T *ivals = (REAL_T *) A->data_[i]->vals_;
        const int innz = A->data_[i]->NNZ_;
        LOG_INFO("annz = %d, j0 = %d, jn = %d\n", innz, icols[0],
                 icols[innz - 1]);
//        for (int ipos = innz - 1; ipos >= 0; ipos--)
        int ipos = innz - 1;
        while (ipos >= nz_t[i])
        {
            const int j = icols[ipos];

            if (j > i)
            {
                int *jcols = A->data_[j]->cols_;
                REAL_T *jvals = (REAL_T *) A->data_[j]->vals_;
                const int jnnz = A->data_[j]->NNZ_;
                const int jstart = nz_t[j];
                int found = 0;
                /* search for symmetric position */
                for (int jpos = jstart; jpos < jnnz; jpos++)
                {
                    const int k = jcols[jpos];
                    if (k == i)
                    {
                        /* symmetric position found so just swap entries */
                        REAL_T tmp = ivals[ipos];
                        ivals[ipos] = jvals[jpos];
                        jvals[jpos] = tmp;
                        /* swap position in row i */
                        TYPED_FUNC(csr_swap_row_entries) (A->data_[i], ipos,
                                                          nz_t[i]);
                        /* swap position in row j */
                        TYPED_FUNC(csr_swap_row_entries) (A->data_[j], jpos,
                                                          nz_t[j]);
                        /* update nonzero count */
                        nz_t[i]++;
                        nz_t[j]++;
                        found = 1;
                        break;
                    }
                }
                if (!found)
                {
                    /* nonsymmetric entry. Insert entry and swap position */
                    TYPED_FUNC(csr_set_row_element_new) (A->data_[j], i,
                                                         &ivals[ipos]);
                    /* swap position in row j */
                    const int nzpos = csr_row_NNZ(A->data_[j]) - 1;
                    TYPED_FUNC(csr_swap_row_entries) (A->data_[j], nzpos,
                                                      nz_t[j]);
                    /* update nonzero count */
                    nz_t[j]++;
                    /* update ipos */
                    ipos--;
                }
            }
            else
            {
                if (j == i)
                {
                    LOG_INFO("ipos = %d, j = %d, nz_t = %d\n", ipos, j,
                             nz_t[i]);
                    /* swap position in row i */
                    TYPED_FUNC(csr_swap_row_entries) (A->data_[i], ipos,
                                                      nz_t[i]);
                    /* update nonzero count */
                    nz_t[i]++;
                }
                else
                {
                    ipos--;
                }
            }
        }
    }
    // nonsymmetric contribution from lower triangular part
    for (int i = 0; i < N; i++)
    {
        int *icols = A->data_[i]->cols_;
        REAL_T *ivals = (REAL_T *) A->data_[i]->vals_;
        const int innz = A->data_[i]->NNZ_;
        LOG_INFO("nz_t[i] = %d, nnz = %d\n", nz_t[i], innz);
        for (int ipos = nz_t[i]; ipos < innz; ipos++)
        {
            const int j = icols[ipos];
            if (j < i)
            {
                // insert entry in row j
                TYPED_FUNC(csr_set_row_element_new) (A->data_[j], i,
                                                     &ivals[ipos]);
                /* swap position in row j */
                const int nzpos = csr_row_NNZ(A->data_[j]) - 1;
                TYPED_FUNC(csr_swap_row_entries) (A->data_[j], nzpos,
                                                  nz_t[j]);
                /* update nonzero count */
                nz_t[j]++;
            }
        }
    }
    /* update nonzero row counts */
#pragma omp parallel for shared(N)
    for (int i = 0; i < N; i++)
    {
        csr_row_NNZ(A->data_[i]) = nz_t[i];
    }
}
