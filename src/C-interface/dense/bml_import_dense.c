#include "../bml_allocate.h"
#include "../bml_logger.h"
#include "../bml_types.h"
#include "bml_allocate_dense.h"
#include "bml_import_dense.h"
#include "bml_types_dense.h"

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
bml_matrix_dense_t *
bml_import_from_dense_dense(
    bml_matrix_precision_t matrix_precision,
    bml_dense_order_t order,
    int N,
    void *A,
    double threshold,
    bml_distribution_mode_t distrib_mode)
{
    bml_matrix_dense_t *A_bml = NULL;

    switch (matrix_precision)
    {
        case single_real:
            A_bml =
                bml_import_from_dense_dense_single_real(order, N, A,
                                                        distrib_mode);
            break;
        case double_real:
            A_bml =
                bml_import_from_dense_dense_double_real(order, N, A,
                                                        distrib_mode);
            break;
        case single_complex:
            A_bml =
                bml_import_from_dense_dense_single_complex(order, N, A,
                                                           distrib_mode);
            break;
        case double_complex:
            A_bml =
                bml_import_from_dense_dense_double_complex(order, N, A,
                                                           distrib_mode);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
    return A_bml;
}
