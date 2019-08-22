#include "../bml_add.h"
#include "../bml_logger.h"
#include "../bml_multiply.h"
#include "../bml_types.h"
#include "bml_add_ellsort.h"
#include "bml_multiply_ellsort.h"
#include "bml_types_ellsort.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/** Matrix multiply.
 *
 * C = alpha * A * B + beta * C
 *
 *  \ingroup multiply_group
 *
 *  \param A Matrix A
 *  \param B Matrix B
 *  \param C Matrix C
 *  \param alpha Scalar factor multiplied by A * B
 *  \param beta Scalar factor multiplied by C
 *  \param threshold Used for sparse multiply
 */
void
bml_multiply_ellsort(
    bml_matrix_ellsort_t * A,
    bml_matrix_ellsort_t * B,
    bml_matrix_ellsort_t * C,
    double alpha,
    double beta,
    double threshold)
{
    switch (A->matrix_precision)
    {
        case single_real:
            bml_multiply_ellsort_single_real(A, B, C, alpha, beta, threshold);
            break;
        case double_real:
            bml_multiply_ellsort_double_real(A, B, C, alpha, beta, threshold);
            break;
        case single_complex:
            bml_multiply_ellsort_single_complex(A, B, C, alpha, beta,
                                                threshold);
            break;
        case double_complex:
            bml_multiply_ellsort_double_complex(A, B, C, alpha, beta,
                                                threshold);
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
 *  \param threshold Used for sparse multiply
 */
void *
bml_multiply_x2_ellsort(
    bml_matrix_ellsort_t * X,
    bml_matrix_ellsort_t * X2,
    double threshold)
{
    switch (X->matrix_precision)
    {
        case single_real:
            return bml_multiply_x2_ellsort_single_real(X, X2, threshold);
            break;
        case double_real:
            return bml_multiply_x2_ellsort_double_real(X, X2, threshold);
            break;
        case single_complex:
            return bml_multiply_x2_ellsort_single_complex(X, X2, threshold);
            break;
        case double_complex:
            return bml_multiply_x2_ellsort_double_complex(X, X2, threshold);
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
 *  \param A Matrix A
 *  \param B Matrix B
 *  \param C Matrix C
 *  \param threshold Used for sparse multiply
 */
void
bml_multiply_AB_ellsort(
    bml_matrix_ellsort_t * A,
    bml_matrix_ellsort_t * B,
    bml_matrix_ellsort_t * C,
    double threshold)
{
    switch (A->matrix_precision)
    {
        case single_real:
            bml_multiply_AB_ellsort_single_real(A, B, C, threshold);
            break;
        case double_real:
            bml_multiply_AB_ellsort_double_real(A, B, C, threshold);
            break;
        case single_complex:
            bml_multiply_AB_ellsort_single_complex(A, B, C, threshold);
            break;
        case double_complex:
            bml_multiply_AB_ellsort_double_complex(A, B, C, threshold);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
}

/** Matrix multiply with threshold adjustment.
 *
 * C = A * B
 *
 *  \ingroup multiply_group
 *
 *  \param A Matrix A
 *  \param B Matrix B
 *  \param C Matrix C
 *  \param threshold Used for sparse multiply
 */
void
bml_multiply_adjust_AB_ellsort(
    bml_matrix_ellsort_t * A,
    bml_matrix_ellsort_t * B,
    bml_matrix_ellsort_t * C,
    double threshold)
{
    switch (A->matrix_precision)
    {
        case single_real:
            bml_multiply_adjust_AB_ellsort_single_real(A, B, C, threshold);
            break;
        case double_real:
            bml_multiply_adjust_AB_ellsort_double_real(A, B, C, threshold);
            break;
        case single_complex:
            bml_multiply_adjust_AB_ellsort_single_complex(A, B, C, threshold);
            break;
        case double_complex:
            bml_multiply_adjust_AB_ellsort_double_complex(A, B, C, threshold);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
}
