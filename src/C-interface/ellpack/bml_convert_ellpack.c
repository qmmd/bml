#include "bml_convert_ellpack.h"
#include "bml_logger.h"

#include <stdlib.h>

bml_matrix_ellpack_t *
bml_convert_ellpack(
    const bml_matrix_t * A,
    const bml_matrix_precision_t matrix_precision,
    const int M,
    const bml_distribution_mode_t distrib_mode)
{
    switch (matrix_precision)
    {
        case single_real:
            return bml_convert_ellpack_single_real(A, matrix_precision, M, distrib_mode);
            break;
        case double_real:
            return bml_convert_ellpack_double_real(A, matrix_precision, M, distrib_mode);
            break;
        case single_complex:
            return bml_convert_ellpack_single_complex(A, matrix_precision, M, distrib_mode);
            break;
        case double_complex:
            return bml_convert_ellpack_double_complex(A, matrix_precision, M, distrib_mode);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
    return NULL;
}
