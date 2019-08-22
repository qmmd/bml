#include "../bml_allocate.h"
#include "../bml_logger.h"
#include "bml_allocate_ellblock.h"
#include "bml_import_ellblock.h"
#include "bml_types_ellblock.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

/** Convert a dense matrix into a bml matrix.
 *
 * \ingroup convert_group
 *
 * \param N The number of rows/columns
 * \param matrix_precision The real precision
 * \param A The dense matrix
 * \param threshold The matrix element magnited threshold
 * \return The bml matrix
 */
bml_matrix_ellblock_t *
bml_import_from_dense_ellblock(
    bml_matrix_precision_t matrix_precision,
    bml_dense_order_t order,
    int N,
    void *A,
    double threshold,
    int M,
    bml_distribution_mode_t distrib_mode)
{
    switch (matrix_precision)
    {
        case single_real:
            return bml_import_from_dense_ellblock_single_real(order, N, A,
                                                              threshold, M,
                                                              distrib_mode);
            break;
        case double_real:
            return bml_import_from_dense_ellblock_double_real(order, N, A,
                                                              threshold, M,
                                                              distrib_mode);
            break;
        case single_complex:
            return bml_import_from_dense_ellblock_single_complex(order, N, A,
                                                                 threshold,
                                                                 M,
                                                                 distrib_mode);
            break;
        case double_complex:
            return bml_import_from_dense_ellblock_double_complex(order, N, A,
                                                                 threshold,
                                                                 M,
                                                                 distrib_mode);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
    return NULL;
}
