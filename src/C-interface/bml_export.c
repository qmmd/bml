#include "bml_export.h"
#include "bml_introspection.h"
#include "bml_logger.h"
#include "dense/bml_export_dense.h"
#include "ellpack/bml_export_ellpack.h"
#include "ellsort/bml_export_ellsort.h"
#include "ellblock/bml_export_ellblock.h"
#include "csr/bml_export_csr.h"

#include <stdlib.h>

/** Export a bml matrix.
 *
 * The returned pointer has to be typecase into the proper real
 * type. If the bml matrix is a single precision matrix, then the
 * following should be used:
 *
 * \code{.c}
 * float *A_dense = bml_export_to_dense(A_bml);
 * \endcode
 *
 * The matrix size can be queried with
 *
 * \code{.c}
 * int N = bml_get_size(A_bml);
 * \endcode
 *
 * \ingroup convert_group_C
 *
 * \param A The bml matrix
 * \param order The matrix element order
 * \return The dense matrix
 */
void *
bml_export_to_dense(
    bml_matrix_t * A,
    bml_dense_order_t order)
{
    switch (bml_get_type(A))
    {
        case dense:
            return bml_export_to_dense_dense(A, order);
        case ellpack:
            return bml_export_to_dense_ellpack(A, order);
        case ellsort:
            return bml_export_to_dense_ellsort(A, order);
        case ellblock:
            return bml_export_to_dense_ellblock(A, order);
        case csr:
            return bml_export_to_dense_csr(A, order);
        case type_uninitialized:
            return NULL;
            break;
        default:
            LOG_ERROR("unknown matrix type\n");
    }
    return NULL;
}
