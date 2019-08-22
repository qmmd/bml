#include "../../macros.h"
#include "../../typed.h"
#include "../bml_allocate.h"
#include "../bml_logger.h"
#include "../bml_types.h"
#include "bml_allocate_ellblock.h"
#include "bml_export_ellblock.h"
#include "bml_types_ellblock.h"

#include <assert.h>
#include <complex.h>
#include <math.h>
#include <stdio.h>
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
    bml_export_to_dense_ellblock) (
    bml_matrix_ellblock_t * A,
    bml_dense_order_t order)
{
    assert(A->N > 0);

    int N = A->N;
    int NB = A->NB;
    int MB = A->MB;
    int *A_nnzb = A->nnzb;
    int *A_indexb = A->indexb;
    int *bsize = A->bsize;
    REAL_T *A_dense = bml_allocate_memory(sizeof(REAL_T) * N * N);
    REAL_T **A_ptr_value = (REAL_T **) A->ptr_value;
    int *offset = malloc(NB * sizeof(int));
    offset[0] = 0;
    for (int ib = 1; ib < NB; ib++)
    {
        offset[ib] = offset[ib - 1] + bsize[ib - 1];
        assert(offset[ib] < N);
    }
    switch (order)
    {
        case dense_row_major:
            for (int ib = 0; ib < NB; ib++)
            {
                for (int jp = 0; jp < A_nnzb[ib]; jp++)
                {
                    int ind = ROWMAJOR(ib, jp, NB, MB);
                    int jb = A_indexb[ind];
                    REAL_T *A_value = A_ptr_value[ind];
                    assert(A_value != NULL);
                    for (int ii = 0; ii < bsize[ib]; ii++)
                        for (int jj = 0; jj < bsize[jb]; jj++)
                        {
                            int i = offset[ib] + ii;
                            int j = offset[jb] + jj;
                            A_dense[ROWMAJOR(i, j, N, N)]
                                =
                                A_value[ROWMAJOR
                                        (ii, jj, bsize[ib], bsize[jb])];
                        }
                }
            }
            break;
        case dense_column_major:
            for (int ib = 0; ib < NB; ib++)
            {
                for (int jp = 0; jp < A_nnzb[ib]; jp++)
                {
                    int ind = ROWMAJOR(ib, jp, NB, MB);
                    int jb = A_indexb[ind];
                    REAL_T *A_value = A_ptr_value[ind];
                    for (int ii = 0; ii < bsize[ib]; ii++)
                        for (int jj = 0; jj < bsize[jb]; jj++)
                        {
                            int i = offset[ib] + ii;
                            int j = offset[jb] + jj;
                            A_dense[COLMAJOR(i, j, N, N)]
                                =
                                A_value[ROWMAJOR
                                        (ii, jj, bsize[ib], bsize[jb])];
                        }
                }
            }
            break;
        default:
            LOG_ERROR("unknown order\n");
            break;
    }
    free(offset);

    return A_dense;
}
