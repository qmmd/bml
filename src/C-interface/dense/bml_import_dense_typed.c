#include "../../macros.h"
#include "../../typed.h"
#include "../bml_allocate.h"
#include "../bml_logger.h"
#include "../bml_types.h"
#include "bml_allocate_dense.h"
#include "bml_import_dense.h"
#include "bml_types_dense.h"

#ifdef BML_USE_MAGMA
#include "magma_v2.h"
#endif

#include <complex.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/** Convert a dense matrix into a bml matrix.
 *
 * \ingroup convert_group
 *
 * \param N The number of rows/columns
 * \param matrix_precision The real precision
 * \param A The dense matrix
 * \return The bml matrix
 */
bml_matrix_dense_t *TYPED_FUNC(
    bml_import_from_dense_dense) (
    bml_dense_order_t order,
    int N,
    void *A,
    bml_distribution_mode_t distrib_mode)
{
    bml_matrix_dimension_t matrix_dimension = {N, N, N};
    bml_matrix_dense_t *A_bml =
        TYPED_FUNC(bml_zero_matrix_dense) (matrix_dimension, distrib_mode);
#ifdef BML_USE_MAGMA
    MAGMA(setmatrix) (N, N, A, N, A_bml->matrix, A_bml->ld, A_bml->queue);
#endif
    switch (order)
    {
        case dense_row_major:
#ifndef BML_USE_MAGMA
            memcpy(A_bml->matrix, A, sizeof(REAL_T) * N * N);
#endif
            break;
        case dense_column_major:
        {
#ifdef BML_USE_MAGMA
            MAGMABLAS(transpose_inplace) (N, A_bml->matrix, A_bml->ld,
                                          A_bml->queue);
#else
            REAL_T *A_ptr = (REAL_T *) A;
            REAL_T *B_ptr = (REAL_T *) A_bml->matrix;
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    B_ptr[ROWMAJOR(i, j, N, N)] = A_ptr[COLMAJOR(i, j, N, N)];
                }
            }
#endif
            break;
        }
        default:
            LOG_ERROR("logic error\n");
            break;
    }
    return A_bml;
}
