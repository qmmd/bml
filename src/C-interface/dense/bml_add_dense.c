#include "../bml_add.h"
#include "../bml_logger.h"
#include "../bml_types.h"
#include "bml_add_dense.h"
#include "bml_types_dense.h"

#include <stdlib.h>
#include <string.h>

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
 */
void
bml_add_dense(
    bml_matrix_dense_t * A,
    bml_matrix_dense_t * B,
    double alpha,
    double beta)
{
    switch (A->matrix_precision)
    {
        case single_real:
            bml_add_dense_single_real(A, B, alpha, beta);
            break;
        case double_real:
            bml_add_dense_double_real(A, B, alpha, beta);
            break;
        case single_complex:
            bml_add_dense_single_complex(A, B, alpha, beta);
            break;
        case double_complex:
            bml_add_dense_double_complex(A, B, alpha, beta);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
}

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
 */
double
bml_add_norm_dense(
    bml_matrix_dense_t * A,
    bml_matrix_dense_t * B,
    double alpha,
    double beta)
{
    double trnorm = 0.0;

    switch (A->matrix_precision)
    {
        case single_real:
            trnorm = bml_add_norm_dense_single_real(A, B, alpha, beta);
            break;
        case double_real:
            trnorm = bml_add_norm_dense_double_real(A, B, alpha, beta);
            break;
        case single_complex:
            trnorm = bml_add_norm_dense_single_complex(A, B, alpha, beta);
            break;
        case double_complex:
            trnorm = bml_add_norm_dense_double_complex(A, B, alpha, beta);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
    return trnorm;
}

/** Matrix addition.
 *
 * \f$ A = A + \beta \mathrm{Id} \f$
 *
 *  \ingroup add_group
 *
 *  \param A Matrix A
 *  \param beta Scalar factor multiplied by I
 */
void
bml_add_identity_dense(
    bml_matrix_dense_t * A,
    double beta)
{
    switch (A->matrix_precision)
    {
        case single_real:
            bml_add_identity_dense_single_real(A, beta);
            break;
        case double_real:
            bml_add_identity_dense_double_real(A, beta);
            break;
        case single_complex:
            bml_add_identity_dense_single_complex(A, beta);
            break;
        case double_complex:
            bml_add_identity_dense_double_complex(A, beta);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
}

/** Matrix addition.
 *
 * \f$ A = \alpha A + \beta \mathrm{Id} \f$
 *
 *  \ingroup add_group
 *
 *  \param A Matrix A
 *  \param alpha Scalar factor multiplied by A
 *  \param beta Scalar factor multiplied by I
 */
void
bml_scale_add_identity_dense(
    bml_matrix_dense_t * A,
    double alpha,
    double beta)
{
    switch (A->matrix_precision)
    {
        case single_real:
            bml_scale_add_identity_dense_single_real(A, alpha, beta);
            break;
        case double_real:
            bml_scale_add_identity_dense_double_real(A, alpha, beta);
            break;
        case single_complex:
            bml_scale_add_identity_dense_single_complex(A, alpha, beta);
            break;
        case double_complex:
            bml_scale_add_identity_dense_double_complex(A, alpha, beta);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
}
