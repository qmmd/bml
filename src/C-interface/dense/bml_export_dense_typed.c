#include "../../macros.h"
#include "../../typed.h"
#include "../bml_allocate.h"
#include "../bml_logger.h"
#include "../bml_types.h"
#include "bml_allocate_dense.h"
#include "bml_export_dense.h"
#include "bml_types_dense.h"

#ifdef BML_USE_MAGMA
#include "magma_v2.h"
#endif

#include <complex.h>
#include <stdlib.h>
#include <string.h>

/** Convert a bml matrix into a dense matrix.
 *
 * \ingroup convert_group
 *
 * \param A The bml matrix
 * \return The dense matrix
 */
void *TYPED_FUNC(
    bml_export_to_dense_dense) (
    bml_matrix_dense_t * A,
    bml_dense_order_t order)
{
    REAL_T *A_dense = bml_allocate_memory(sizeof(REAL_T) * A->N * A->N);

    switch (order)
    {
        case dense_row_major:
#ifdef BML_USE_MAGMA
            MAGMA(getmatrix) (A->N, A->N,
                              A->matrix, A->ld,
                              (MAGMA_T *) A_dense, A->N, A->queue);
#else
            memcpy(A_dense, A->matrix, sizeof(REAL_T) * A->N * A->N);
#endif
            break;
        case dense_column_major:
        {
#ifdef BML_USE_MAGMA
            MAGMABLAS(transpose_inplace) (A->N, A->matrix, A->ld, A->queue);
            MAGMA(getmatrix) (A->N, A->N,
                              A->matrix, A->ld,
                              (MAGMA_T *) A_dense, A->N, A->queue);
            MAGMABLAS(transpose_inplace) (A->N, A->matrix, A->ld, A->queue);
#else
            REAL_T *B_ptr = (REAL_T *) A->matrix;
            for (int i = 0; i < A->N; i++)
            {
                for (int j = 0; j < A->N; j++)
                {
                    A_dense[COLMAJOR(i, j, A->N, A->N)] =
                        B_ptr[ROWMAJOR(i, j, A->N, A->N)];
                }
            }
#endif
            break;
        }
        default:
            LOG_ERROR("logic error\n");
            break;
    }
    return A_dense;
}
