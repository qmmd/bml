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
#include "bml_scale_csr.h"
#include "bml_introspection_csr.h"
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
    int tsize = bml_get_bandwidth_csr(A);
#pragma omp parallel default(none) \
    shared(N, tsize, A, B) \
    shared(alpha, beta, threshold)
    {
        /* create hash table */
        csr_row_index_hash_t *table = csr_noinit_table(tsize);
#pragma omp for
        for (int i = 0; i < N; i++)
        {
            int *acols = A->data_[i]->cols_;
            REAL_T *avals = (REAL_T *) A->data_[i]->vals_;
            const int annz = A->data_[i]->NNZ_;

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
            //reset table
            csr_reset_table(table);
        }
        // delete table
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
    int N = A->N_;

#pragma omp parallel for                  \
    shared(N)
    for (int i = 0; i < N; i++)
    {
        int *acols = A->data_[i]->cols_;
        REAL_T *avals = (REAL_T *) A->data_[i]->vals_;
        const int annz = A->data_[i]->NNZ_;

        int diag = -1;

        // find position of diagonal entry
        for (int pos = 0; pos < annz; pos++)
        {
            if (acols[pos] == i)
            {
                diag = pos;
                break;
            }
        }

        if (beta > (double) 0.0 || beta < (double) 0.0)
        {
            // if diagonal entry does not exist, insert, else add
            REAL_T val = (REAL_T) beta;
            if (diag == -1)
            {
                TYPED_FUNC(csr_set_row_element_new) (A->data_[i], i, &val);
            }
            else
            {
                avals[diag] += val;
            }

        }
    }
    /* apply thresholding */
    TYPED_FUNC(bml_threshold_csr) (A, threshold);
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
    // scale then update diagonal
    TYPED_FUNC(bml_scale_inplace_csr) (&alpha, A);

    TYPED_FUNC(bml_add_identity_csr) (A, beta, threshold);
}
