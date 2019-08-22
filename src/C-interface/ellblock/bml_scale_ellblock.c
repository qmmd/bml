#include "../bml_logger.h"
#include "../bml_scale.h"
#include "../bml_types.h"
#include "bml_scale_ellblock.h"
#include "bml_types_ellblock.h"

#include <stdlib.h>
#include <string.h>

/** Scale an ellblock matrix - result is a new matrix.
 *
 *  \ingroup scale_group
 *
 *  \param A The matrix to be scaled
 *  \return A scale version of matrix A.
 */
bml_matrix_ellblock_t *
bml_scale_ellblock_new(
    void *scale_factor,
    bml_matrix_ellblock_t * A)
{
    bml_matrix_ellblock_t *B = NULL;

    switch (A->matrix_precision)
    {
        case single_real:
            B = bml_scale_ellblock_new_single_real(scale_factor, A);
            break;
        case double_real:
            B = bml_scale_ellblock_new_double_real(scale_factor, A);
            break;
        case single_complex:
            B = bml_scale_ellblock_new_single_complex(scale_factor, A);
            break;
        case double_complex:
            B = bml_scale_ellblock_new_double_complex(scale_factor, A);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
    return B;
}

/** Scale an ellblock matrix.
 *
 *  \ingroup scale_group
 *
 *  \param A The matrix to be scaled
 *  \param B Scaled version of matrix A
 */
void
bml_scale_ellblock(
    void *scale_factor,
    bml_matrix_ellblock_t * A,
    bml_matrix_ellblock_t * B)
{
    switch (A->matrix_precision)
    {
        case single_real:
            bml_scale_ellblock_single_real(scale_factor, A, B);
            break;
        case double_real:
            bml_scale_ellblock_double_real(scale_factor, A, B);
            break;
        case single_complex:
            bml_scale_ellblock_single_complex(scale_factor, A, B);
            break;
        case double_complex:
            bml_scale_ellblock_double_complex(scale_factor, A, B);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
}

void
bml_scale_inplace_ellblock(
    void *scale_factor,
    bml_matrix_ellblock_t * A)
{
    switch (A->matrix_precision)
    {
        case single_real:
            bml_scale_inplace_ellblock_single_real(scale_factor, A);
            break;
        case double_real:
            bml_scale_inplace_ellblock_double_real(scale_factor, A);
            break;
        case single_complex:
            bml_scale_inplace_ellblock_single_complex(scale_factor, A);
            break;
        case double_complex:
            bml_scale_inplace_ellblock_double_complex(scale_factor, A);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
}
