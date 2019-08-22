#include "../bml_inverse.h"
#include "../bml_logger.h"
#include "../bml_types.h"
#include "bml_inverse_dense.h"
#include "bml_types_dense.h"

#include <stdlib.h>
#include <string.h>

/** Matrix inverse.
 *
 * \ingroup inverse_group
 *
 * \param A Matrix A
 * \return Matrix inverse of A
 */
bml_matrix_dense_t *
bml_inverse_dense(
    bml_matrix_dense_t * A)
{

    bml_matrix_dense_t *B = NULL;

    switch (A->matrix_precision)
    {
        case single_real:
            B = bml_inverse_dense_single_real(A);
            break;
        case double_real:
            B = bml_inverse_dense_double_real(A);
            break;
        case single_complex:
            B = bml_inverse_dense_single_complex(A);
            break;
        case double_complex:
            B = bml_inverse_dense_double_complex(A);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }

    return B;

}

/** Matrix inverse inplace.
 *
 * \ingroup inverse_group
 *
 * \param A Matrix A
 * \param B Matrix inverse of A
 */
void
bml_inverse_inplace_dense(
    bml_matrix_dense_t * A)
{
    switch (A->matrix_precision)
    {
        case single_real:
            bml_inverse_inplace_dense_single_real(A);
            break;
        case double_real:
            bml_inverse_inplace_dense_double_real(A);
            break;
        case single_complex:
            bml_inverse_inplace_dense_single_complex(A);
            break;
        case double_complex:
            bml_inverse_inplace_dense_double_complex(A);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
}
