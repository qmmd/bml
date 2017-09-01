#include "bml_convert.h"
#include "bml_logger.h"
#include "dense/bml_convert_dense.h"

#include <stdlib.h>

/** Convert a bml matrix to another type.
 *
 * \f$ A \rightarrow B \f$
 *
 * \param A The input matrix.
 * \return The converted matrix \f$ B \f$.
 */
bml_matrix_t *
bml_convert(
    const bml_matrix_t * A,
    const bml_matrix_type_t matrix_type,
    const bml_matrix_precision_t matrix_precision,
    const int M,
    const bml_distribution_mode_t distrib_mode)
{
    switch (matrix_type)
    {
        case dense:
            return bml_convert_dense(A, matrix_precision, distrib_mode);
            break;
        default:
            LOG_ERROR("unknown matrix type\n");
            break;
    }
    return NULL;
}
