#include "../bml_allocate.h"
#include "../bml_logger.h"
#include "../bml_scale.h"
#include "../bml_types.h"
#include "bml_allocate_dense.h"
#include "bml_copy_dense.h"
#include "bml_scale_dense.h"
#include "bml_types_dense.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

/** Scale a dense matrix - result in new matrix.
 *
 * \ingroup scale_group
 *
 * \param A The matrix to be scaled
 * \return A scaled version of matrix A.
 */
bml_matrix_dense_t *
bml_scale_dense_new(
    void *scale_factor,
    bml_matrix_dense_t * A)
{
    bml_matrix_dense_t *B = NULL;

    assert(A != NULL);

    switch (A->matrix_precision)
    {
        case single_real:
            B = bml_scale_dense_new_single_real(scale_factor, A);
            break;
        case double_real:
            B = bml_scale_dense_new_double_real(scale_factor, A);
            break;
        case single_complex:
            B = bml_scale_dense_new_single_complex(scale_factor, A);
            break;
        case double_complex:
            B = bml_scale_dense_new_double_complex(scale_factor, A);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
    return B;
}

/** Scale a dense matrix.
 *
 *  \ingroup scale_group
 *
 *  \param A The matrix to be scaled
 *  \param B Scaled version of matrix A
 */
void
bml_scale_dense(
    void *scale_factor,
    bml_matrix_dense_t * A,
    bml_matrix_dense_t * B)
{
    assert(A != NULL);

    switch (A->matrix_precision)
    {
        case single_real:
            bml_scale_dense_single_real(scale_factor, A, B);
            break;
        case double_real:
            bml_scale_dense_double_real(scale_factor, A, B);
            break;
        case single_complex:
            bml_scale_dense_single_complex(scale_factor, A, B);
            break;
        case double_complex:
            bml_scale_dense_double_complex(scale_factor, A, B);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
}

void
bml_scale_inplace_dense(
    void *scale_factor,
    bml_matrix_dense_t * A)
{
    switch (A->matrix_precision)
    {
        case single_real:
            bml_scale_inplace_dense_single_real(scale_factor, A);
            break;
        case double_real:
            bml_scale_inplace_dense_double_real(scale_factor, A);
            break;
        case single_complex:
            bml_scale_inplace_dense_single_complex(scale_factor, A);
            break;
        case double_complex:
            bml_scale_inplace_dense_double_complex(scale_factor, A);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
}
