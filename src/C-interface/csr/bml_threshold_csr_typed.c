#include "../../macros.h"
#include "../../typed.h"
#include "../bml_allocate.h"
#include "../bml_parallel.h"
#include "../bml_threshold.h"
#include "../bml_types.h"
#include "bml_allocate_csr.h"
#include "bml_threshold_csr.h"
#include "bml_types_csr.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/** Threshold a matrix.
 *
 *  \ingroup threshold_group
 *
 *  \param A The matrix to be thresholded
 *  \param threshold Threshold value
 *  \return the thresholded A
 */
bml_matrix_csr_t
    * TYPED_FUNC(bml_threshold_new_csr) (bml_matrix_csr_t * A,
                                         double threshold)
{
    int N = A->N_;
    int M = A->NZMAX_;

    bml_matrix_csr_t *B =
        TYPED_FUNC(bml_zero_matrix_csr) (N, M, A->distribution_mode);

#pragma omp parallel for               \
    shared(N)
    for (int i = 0; i < N; i++)
    {
        int *bcols = B->data_[i]->cols_;
        REAL_T *bvals = (REAL_T *) B->data_[i]->vals_;
        const int bnnz = B->data_[i]->NNZ_;

        int *cols = A->data_[i]->cols_;
        REAL_T *vals = (REAL_T *) A->data_[i]->vals_;
        const int annz = A->data_[i]->NNZ_;
        for (int pos = 0; pos < annz; pos++)
        {
            if (is_above_threshold(vals[pos], threshold))
            {
//                bml_set_element_new_csr(B, i, cols[pos], &vals[pos]);
                bvals[bnnz] = vals[pos];
                bcols[bnnz] = cols[pos];
                /* update nnz of row i of B */
                B->data_[i]->NNZ_++;
            }
        }
    }

    return B;
}

/** Threshold a matrix in place.
 *
 *  \ingroup threshold_group
 *
 *  \param A The matrix to be thresholded
 *  \param threshold Threshold value
 *  \return the thresholded A
 */
void TYPED_FUNC(
    bml_threshold_csr) (
    bml_matrix_csr_t * A,
    double threshold)
{
    int N = A->N_;

    int rlen;
#pragma omp parallel for               \
    private(rlen) \
    shared(N)
    for (int i = 0; i < N; i++)
    {
        rlen = 0;
        int *cols = A->data_[i]->cols_;
        REAL_T *vals = (REAL_T *) A->data_[i]->vals_;
        const int annz = A->data_[i]->NNZ_;
        for (int pos = 0; pos < annz; pos++)
        {
            const int j = cols[pos];
            if (is_above_threshold(vals[pos], threshold))
            {
                if (rlen < j)
                {
                    /* move row entry */
                    cols[rlen] = cols[pos];
                    vals[rlen] = vals[pos];
                }
                rlen++;
            }
        }
        csr_row_NNZ(A->data_[i]) = rlen;
    }
}
