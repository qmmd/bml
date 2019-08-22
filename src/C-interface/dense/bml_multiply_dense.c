#include "../bml_logger.h"
#include "../bml_multiply.h"
#include "../bml_types.h"
#include "bml_multiply_dense.h"
#include "bml_types_dense.h"

#include <stdlib.h>
#include <string.h>

/** Matrix multiply.
 *
 * \f$ C \leftarrow \alpha \, A \, B + \beta C \f$
 *
 * \ingroup multiply_group
 *
 * \param A Matrix A
 * \param B Matrix B
 * \param C Matrix C
 * \param alpha Scalar factor multiplied by A * B
 * \param beta Scalar factor multiplied by C
 */
void
bml_multiply_dense(
    bml_matrix_dense_t * A,
    bml_matrix_dense_t * B,
    bml_matrix_dense_t * C,
    double alpha,
    double beta)
{
    switch (A->matrix_precision)
    {
        case single_real:
            bml_multiply_dense_single_real(A, B, C, alpha, beta);
            break;
        case double_real:
            bml_multiply_dense_double_real(A, B, C, alpha, beta);
            break;
        case single_complex:
            bml_multiply_dense_single_complex(A, B, C, alpha, beta);
            break;
        case double_complex:
            bml_multiply_dense_double_complex(A, B, C, alpha, beta);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
}

/** Matrix multiply.
 *
 * X2 = X * X
 *
 *  \ingroup multiply_group
 *
 *  \param X Matrix X
 *  \param X2 Matrix X2
 */
void *
bml_multiply_x2_dense(
    bml_matrix_dense_t * X,
    bml_matrix_dense_t * X2)
{
    switch (X->matrix_precision)
    {
        case single_real:
            return bml_multiply_x2_dense_single_real(X, X2);
            break;
        case double_real:
            return bml_multiply_x2_dense_double_real(X, X2);
            break;
        case single_complex:
            return bml_multiply_x2_dense_single_complex(X, X2);
            break;
        case double_complex:
            return bml_multiply_x2_dense_double_complex(X, X2);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
    return NULL;
}

/** Matrix multiply.
 *
 * C = A * B
 *
 *  \ingroup multiply_group
 *
 *   \param A Matrix A
 *   \param B Matrix B
 *   \param C Matrix C
 */
void
bml_multiply_AB_dense(
    bml_matrix_dense_t * A,
    bml_matrix_dense_t * B,
    bml_matrix_dense_t * C)
{
    switch (A->matrix_precision)
    {
        case single_real:
            bml_multiply_AB_dense_single_real(A, B, C);
            break;
        case double_real:
            bml_multiply_AB_dense_double_real(A, B, C);
            break;
        case single_complex:
            bml_multiply_AB_dense_single_complex(A, B, C);
            break;
        case double_complex:
            bml_multiply_AB_dense_double_complex(A, B, C);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
}

/** Matrix multiply.
 *
 * Note: This routine is provided for completeness.
 *
 * C = A * B
 *
 *  \ingroup multiply_group
 *
 *   \param A Matrix A
 *   \param B Matrix B
 *   \param C Matrix C
 */
void
bml_multiply_adjust_AB_dense(
    bml_matrix_dense_t * A,
    bml_matrix_dense_t * B,
    bml_matrix_dense_t * C)
{
    switch (A->matrix_precision)
    {
        case single_real:
            bml_multiply_AB_dense_single_real(A, B, C);
            break;
        case double_real:
            bml_multiply_AB_dense_double_real(A, B, C);
            break;
        case single_complex:
            bml_multiply_AB_dense_single_complex(A, B, C);
            break;
        case double_complex:
            bml_multiply_AB_dense_double_complex(A, B, C);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
}
