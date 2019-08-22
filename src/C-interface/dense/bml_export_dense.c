#include "../bml_allocate.h"
#include "../bml_logger.h"
#include "../bml_types.h"
#include "bml_allocate_dense.h"
#include "bml_export_dense.h"
#include "bml_types_dense.h"

#ifdef BML_USE_MAGMA
#include "magma_v2.h"
#endif

#include <stdlib.h>
#include <string.h>

/** Convert a bml matrix into a dense matrix.
 *
 * \ingroup convert_group
 *
 * \param A The bml matrix
 * \return The dense matrix
 */
void *
bml_export_to_dense_dense(
    bml_matrix_dense_t * A,
    bml_dense_order_t order)
{
    switch (A->matrix_precision)
    {
        case single_real:
            return bml_export_to_dense_dense_single_real(A, order);
            break;
        case double_real:
            return bml_export_to_dense_dense_double_real(A, order);
            break;
        case single_complex:
            return bml_export_to_dense_dense_single_complex(A, order);
            break;
        case double_complex:
            return bml_export_to_dense_dense_double_complex(A, order);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
    return NULL;
}
