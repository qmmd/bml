#include "../bml_logger.h"
#include "../bml_transpose_triangle.h"
#include "../bml_types.h"
#include "bml_transpose_triangle_dense.h"
#include "bml_types_dense.h"

#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif


/** Transposes a triangle of a matrix in place.
 *
 *  \ingroup transpose_triangle_group
 *
 *  \param A  The matrix for which the triangle should be transposed
 *  \param triangle  Which triangle to transpose ('u': upper, 'l': lower)
 */
void
bml_transpose_triangle_dense(
    bml_matrix_dense_t * A,
    char triangle)
{
    switch (A->matrix_precision)
    {
        case single_real:
            bml_transpose_triangle_dense_single_real(A, triangle);
            break;
        case double_real:
            bml_transpose_triangle_dense_double_real(A, triangle);
            break;
        case single_complex:
            bml_transpose_triangle_dense_single_complex(A, triangle);
            break;
        case double_complex:
            bml_transpose_triangle_dense_double_complex(A, triangle);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
}
