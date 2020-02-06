#include "../../macros.h"
#include "../../typed.h"
#include "../bml_allocate.h"
#include "../bml_logger.h"
#include "../bml_types.h"
#include "bml_allocate_ellpack.h"
#include "bml_export_ellpack.h"
#include "bml_types_ellpack.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/** Convert a bml matrix into a dense matrix.
 *
 * \ingroup convert_group
 *
 * \param A The bml matrix
 * \return The dense matrix
 */
void *TYPED_FUNC(
    bml_export_to_dense_ellpack) (
    bml_matrix_ellpack_t * A,
    bml_dense_order_t order)
{
    int N = A->N;
    int M = A->M;
    int *A_nnz = A->nnz;
    int *A_index = A->index;
    REAL_T *A_dense = bml_allocate_memory(sizeof(REAL_T) * A->N * A->N);

#ifdef __INTEL_COMPILER
#pragma omp parallel for simd
    for (int ii = 0; ii < (A->N * A->N); ii++)
    {
        __assume_aligned(A_dense, 64);
        A_dense[ii] = 0.0;
    }
#endif

    REAL_T *A_value = A->value;

    switch (order)
    {
        case dense_row_major:
#pragma omp parallel for shared(N, M, A_nnz, A_index, A_value, A_dense)
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < A_nnz[i]; j++)
                {
                    A_dense[ROWMAJOR
                            (i, A_index[ROWMAJOR(i, j, N, M)], N,
                             N)] = A_value[ROWMAJOR(i, j, N, M)];
                }
            }
            break;
        case dense_column_major:
#pragma omp parallel for shared(N, M, A_nnz, A_index, A_value, A_dense)
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < A_nnz[i]; j++)
                {
                    A_dense[COLMAJOR
                            (i, A_index[ROWMAJOR(i, j, N, M)], N,
                             N)] = A_value[ROWMAJOR(i, j, N, M)];
                }
            }
            break;
        default:
            LOG_ERROR("unknown order\n");
            break;
    }
    return A_dense;
}
