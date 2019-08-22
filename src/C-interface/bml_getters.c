#include "bml_introspection.h"
#include "bml_logger.h"
#include "bml_getters.h"
#include "dense/bml_getters_dense.h"
#include "ellpack/bml_getters_ellpack.h"
#include "ellsort/bml_getters_ellsort.h"
#include "ellblock/bml_getters_ellblock.h"

#include <stdio.h>

/** Return a single matrix element.
 *
 * \param i The row index
 * \param j The column index
 * \param A The bml matrix
 * \return The matrix element
 */
void *
bml_get(
    bml_matrix_t * A,
    int i,
    int j)
{
    switch (bml_get_type(A))
    {
        case dense:
            return bml_get_dense(A, i, j);
            break;
        case ellpack:
            return bml_get_ellpack(A, i, j);
            break;
        case ellsort:
            return bml_get_ellsort(A, i, j);
            break;
        case ellblock:
            return bml_get_ellblock(A, i, j);
            break;
        default:
            LOG_ERROR("unknown matrix type\n");
            break;
    }
    return NULL;
}

/** Get a whole row.
 *
 * @param A The matrix.
 * @param i The row index.
 * @return An array (needs to be cast into the appropriate type).
 */
void *
bml_get_row(
    bml_matrix_t * A,
    int i)
{
    switch (bml_get_type(A))
    {
        case dense:
            return bml_get_row_dense(A, i);
            break;
        case ellpack:
            return bml_get_row_ellpack(A, i);
            break;
        case ellsort:
            return bml_get_row_ellsort(A, i);
            break;
        case ellblock:
            return bml_get_row_ellblock(A, i);
            break;
        default:
            LOG_ERROR("unknown matrix type in bml_get_row\n");
            break;
    }
    return NULL;
}

/** Get the diagonal.
 *
 * @param A The matrix.
 * @return The diagonal (an array)
 */
void *
bml_get_diagonal(
    bml_matrix_t * A)
{
    switch (bml_get_type(A))
    {
        case dense:
            return bml_get_diagonal_dense(A);
            break;
        case ellpack:
            return bml_get_diagonal_ellpack(A);
            break;
        case ellsort:
            return bml_get_diagonal_ellsort(A);
            break;
        case ellblock:
            return bml_get_diagonal_ellblock(A);
            break;
        default:
            LOG_ERROR("unknown matrix type in bml_get_diagonal\n");
            break;
    }
    return NULL;
}
