#include "bml_adjungate_triangle.h"
#include "bml_introspection.h"
#include "bml_logger.h"
#include "dense/bml_adjungate_triangle_dense.h"
#include "ellpack/bml_adjungate_triangle_ellpack.h"

#include <stdlib.h>
#include <stdio.h>

/** Adjungates (conjugate transpose) a triangle of a matrix in place.
 *
 *  \ingroup adjungate_triangle_group
 *
 *  \param A  The matrix for which the triangle should be adjungated
 *  \param triangle  Which triangle to adjungate ('u': upper, 'l': lower)
 */
void
bml_adjungate_triangle(
    bml_matrix_t * A,
    char *triangle)
{
    switch (bml_get_type(A))
    {
        case dense:
            bml_adjungate_triangle_dense(A, triangle);
            break;
        case ellpack:
            bml_adjungate_triangle_ellpack(A, triangle);
            break;
        default:
            LOG_ERROR("unknown matrix type in bml_adjungate_triangle\n");
            break;
    }
}
