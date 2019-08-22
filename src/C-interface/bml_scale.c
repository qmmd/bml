#include "bml_introspection.h"
#include "bml_logger.h"
#include "bml_scale.h"
#include "dense/bml_scale_dense.h"
#include "ellpack/bml_scale_ellpack.h"
#include "ellsort/bml_scale_ellsort.h"
#include "ellblock/bml_scale_ellblock.h"

#include <stdlib.h>

/** Scale a matrix - resulting matrix is new.
 *
 * \ingroup scale_group_C
 *
 * \param scale_factor Scale factor for A
 * \param A Matrix to scale
 * \return A Scaled Copy of A
 */
bml_matrix_t *
bml_scale_new(
    void *scale_factor,
    bml_matrix_t * A)
{
    bml_matrix_t *B = NULL;

    switch (bml_get_type(A))
    {
        case dense:
            B = bml_scale_dense_new(scale_factor, A);
            break;
        case ellpack:
            B = bml_scale_ellpack_new(scale_factor, A);
            break;
        case ellsort:
            B = bml_scale_ellsort_new(scale_factor, A);
            break;
        case ellblock:
            B = bml_scale_ellblock_new(scale_factor, A);
            break;
        default:
            LOG_ERROR("unknown matrix type\n");
            break;
    }
    return B;
}

/** Scale a matrix - resulting matrix exists.
 *
 * \ingroup scale_group_C
 *
 * \param scale_factor Scale factor for A
 * \param A Matrix to scale
 * \param B Scaled Matrix
 */
void
bml_scale(
    void *scale_factor,
    bml_matrix_t * A,
    bml_matrix_t * B)
{
    switch (bml_get_type(A))
    {
        case dense:
            bml_scale_dense(scale_factor, A, B);
            break;
        case ellpack:
            bml_scale_ellpack(scale_factor, A, B);
            break;
        case ellsort:
            bml_scale_ellsort(scale_factor, A, B);
            break;
        case ellblock:
            bml_scale_ellblock(scale_factor, A, B);
            break;
        default:
            LOG_ERROR("unknown matrix type\n");
            break;
    }
}

/** Scale a matrix in place, i.e. the matrix is overwritten.
 *
 * \ingroup scale_group_C
 *
 * \param scale_factor Scale factor for A
 * \param A [inout] Matrix to scale
 */
void
bml_scale_inplace(
    void *scale_factor,
    bml_matrix_t * A)
{
    switch (bml_get_type(A))
    {
        case dense:
            bml_scale_inplace_dense(scale_factor, A);
            break;
        case ellpack:
            bml_scale_inplace_ellpack(scale_factor, A);
            break;
        case ellsort:
            bml_scale_inplace_ellsort(scale_factor, A);
            break;
        case ellblock:
            bml_scale_inplace_ellblock(scale_factor, A);
            break;
        default:
            LOG_ERROR("unknown matrix type\n");
            break;
    }
}
