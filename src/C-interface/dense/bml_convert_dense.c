#include "bml_convert_dense.h"
#include "bml_logger.h"

#include <stdlib.h>

bml_matrix_dense_t *
bml_convert_dense(
    const bml_matrix_t * A,
    const bml_matrix_precision_t matrix_precision,
    const bml_distribution_mode_t distrib_mode)
{
    switch (matrix_precision)
    {
        case single_real:
            return bml_convert_dense_single_real(A, distrib_mode);
            break;
        case double_real:
            return bml_convert_dense_double_real(A, distrib_mode);
            break;
        case single_complex:
            return bml_convert_dense_single_complex(A, distrib_mode);
            break;
        case double_complex:
            return bml_convert_dense_double_complex(A, distrib_mode);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
    return NULL;
}
