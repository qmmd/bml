#include "../../macros.h"
#include "../../typed.h"
#include "../bml_add.h"
#include "../bml_allocate.h"
#include "../bml_parallel.h"
#include "../bml_types.h"
#include "bml_add_csr.h"
#include "bml_allocate_csr.h"
#include "bml_types_csr.h"
#include "bml_setters_csr.h"
#include "bml_threshold_csr.h"
#include "../bml_logger.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/** Matrix addition.
 *
 * \f$ A = \alpha A + \beta B \f$
 *
 * \ingroup add_group
 *
 * \param A Matrix A
 * \param B Matrix B
 * \param alpha Scalar factor multiplied by A
 * \param beta Scalar factor multiplied by B
 * \param threshold Threshold for matrix addition
 */
void TYPED_FUNC(
    bml_add_csr) (
    bml_matrix_csr_t * A,
    bml_matrix_csr_t * B,
    double alpha,
    double beta,
    double threshold)
{
    int N = A->N_;

    for (int i = 0; i < N; i++)
    {
        int *acols = A->data_[i]->cols_;
        REAL_T *avals = (REAL_T *) A->data_[i]->vals_;
        const int annz = A->data_[i]->NNZ_;

        /* create hash table */
        csr_row_index_hash_t *table = csr_noinit_table(annz);
        for (int pos = 0; pos < annz; pos++)
        {
            avals[pos] *= alpha;
            csr_table_insert(table, acols[pos]);
        }
        int *bcols = B->data_[i]->cols_;
        REAL_T *bvals = (REAL_T *) B->data_[i]->vals_;
        const int bnnz = B->data_[i]->NNZ_;
        for (int pos = 0; pos < bnnz; pos++)
        {
            int *idx = (int *) csr_table_lookup(table, bcols[pos]);
            REAL_T val = beta * bvals[pos];
            if (idx)
            {
                avals[*idx] += val;
            }
            else
            {
                TYPED_FUNC(csr_set_row_element_new) (A->data_[i],
                                                     bcols[pos], &val);
            }
        }
        //clear table
        csr_deallocate_table(table);
    }
    /* apply thresholding */
    TYPED_FUNC(bml_threshold_csr) (A, threshold);
}

/******** Not sure why this function is needed or why norms are being computed here -DOK******/
/** Matrix addition.
 *
 * \f$ A = \alpha A + \beta B \f$
 *
 * \ingroup add_group
 *
 * \param A Matrix A
 * \param B Matrix B
 * \param alpha Scalar factor multiplied by A
 * \param beta Scalar factor multiplied by B
 * \param threshold Threshold for matrix addition
 */
double TYPED_FUNC(
    bml_add_norm_csr) (
    bml_matrix_csr_t * A,
    bml_matrix_csr_t * B,
    double alpha,
    double beta,
    double threshold)
{
    LOG_ERROR("bml_add_norm_csr:  Not implemented");
    return 0.;
}

/** Matrix addition.
 *
 *  A = A + beta * I
 *
 *  \ingroup add_group
 *
 *  \param A Matrix A
 *  \param beta Scalar factor multiplied by I
 *  \param threshold Threshold for matrix addition
 */
void TYPED_FUNC(
    bml_add_identity_csr) (
    bml_matrix_csr_t * A,
    double beta,
    double threshold)
{
    REAL_T alpha = (REAL_T) 1.0;

    bml_matrix_csr_t *Id =
        TYPED_FUNC(bml_identity_matrix_csr) (A->N_, A->NZMAX_,
                                             A->distribution_mode);

    TYPED_FUNC(bml_add_csr) (A, Id, alpha, beta, threshold);

    bml_deallocate_csr(Id);
}

/** Matrix addition.
 *
 *  A = alpha * A + beta * I
 *
 *  \ingroup add_group
 *
 *  \param A Matrix A
 *  \param alpha Scalar factor multiplied by A
 *  \param beta Scalar factor multiplied by I
 *  \param threshold Threshold for matrix addition
 */
void TYPED_FUNC(
    bml_scale_add_identity_csr) (
    bml_matrix_csr_t * A,
    double alpha,
    double beta,
    double threshold)
{
    bml_matrix_csr_t *Id =
        TYPED_FUNC(bml_identity_matrix_csr) (A->N_, A->NZMAX_,
                                             A->distribution_mode);

    TYPED_FUNC(bml_add_csr) (A, Id, alpha, beta, threshold);

    bml_deallocate_csr(Id);
}
