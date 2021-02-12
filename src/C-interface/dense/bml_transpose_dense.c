#include "../bml_logger.h"
#include "../bml_transpose.h"
#include "../bml_types.h"
#include "bml_transpose_dense.h"
#include "bml_types_dense.h"

#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/** Threshold a matrix.
 *
 *  \ingroup transpose_group
 *
 *  \param A The matrix to be transposed
 *  \return The transposed A
 */
bml_matrix_dense_t *
bml_transpose_new_dense(
    bml_matrix_dense_t * A)
{
    switch (A->matrix_precision)
    {
        case single_real:
            return bml_transpose_new_dense_single_real(A);
            break;
        case double_real:
            return bml_transpose_new_dense_double_real(A);
            break;
        case single_complex:
            return bml_transpose_new_dense_single_complex(A);
            break;
        case double_complex:
            return bml_transpose_new_dense_double_complex(A);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
    return NULL;
}

/** Transpose a matrix in place.
 *
 *  \ingroup transpose_group
 *
 *  \param A The matrix to be transposed
 *  \return The transposed A
 */
void
bml_transpose_dense(
    bml_matrix_dense_t * A)
{
    switch (A->matrix_precision)
    {
        case single_real:
            bml_transpose_dense_single_real(A);
            break;
        case double_real:
            bml_transpose_dense_double_real(A);
            break;
        case single_complex:
            bml_transpose_dense_single_complex(A);
            break;
        case double_complex:
            bml_transpose_dense_double_complex(A);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
}

/** Complex conjugate a matrix in place.
 *
 * \ingroup transpose_group
 *
 * \param A The matrix to be complex conjugated
 * \return Complex conjugate of A
 */
void
bml_complex_conjugate_dense(
    bml_matrix_dense_t * A)
{
    switch (A->matrix_precision)
    {
        case single_real:
            bml_complex_conjugate_dense_single_real(A);
            break;
        case double_real:
            bml_complex_conjugate_dense_double_real(A);
            break;
        case single_complex:
            bml_complex_conjugate_dense_single_complex(A);
            break;
        case double_complex:
            bml_complex_conjugate_dense_double_complex(A);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
}
