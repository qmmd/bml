#include "../bml_logger.h"
#include "../bml_types.h"
#include "../bml_utilities.h"
#include "../dense/bml_types_dense.h"
#include "bml_diagonalize_ellpack.h"
#include "bml_types_ellpack.h"

#include <string.h>

/** \page diagonalize
 *
 */

void
bml_diagonalize_ellpack(
    const bml_matrix_ellpack_t * A,
    void *eigenvalues,
    bml_matrix_t * eigenvectors)
{
    switch (A->matrix_precision)
    {
        case single_real:
            bml_diagonalize_ellpack_single_real(A, eigenvalues, eigenvectors);
            break;
        case double_real:
            bml_diagonalize_ellpack_double_real(A, eigenvalues, eigenvectors);
            break;
        case single_complex:
            bml_diagonalize_ellpack_single_complex(A, eigenvalues,
                                                   eigenvectors);
            break;
        case double_complex:
            bml_diagonalize_ellpack_double_complex(A, eigenvalues,
                                                   eigenvectors);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
}
