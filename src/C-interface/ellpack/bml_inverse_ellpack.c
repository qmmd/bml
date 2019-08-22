#include "../bml_logger.h"
#include "../bml_types.h"
#include "bml_inverse_ellpack.h"
#include "bml_types_ellpack.h"
#include "../bml_utilities.h"

#include <string.h>

/** \page inverse
 *
 */

bml_matrix_ellpack_t *
bml_inverse_ellpack(
    bml_matrix_ellpack_t * A)
{

    bml_matrix_ellpack_t *B = NULL;

    switch (A->matrix_precision)
    {
        case single_real:
            B = bml_inverse_ellpack_single_real(A);
            break;
        case double_real:
            B = bml_inverse_ellpack_double_real(A);
            break;
        case single_complex:
            B = bml_inverse_ellpack_single_complex(A);
            break;
        case double_complex:
            B = bml_inverse_ellpack_double_complex(A);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }

    return B;

}
