#include "../../macros.h"
#include "../../typed.h"
#include "../bml_allocate.h"
#include "../bml_logger.h"
#include "../bml_types.h"
#include "bml_allocate_ellblock.h"
#include "bml_import_ellblock.h"
#include "bml_types_ellblock.h"
#include "bml_utilities_ellblock.h"

#include <assert.h>
#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/** Convert a dense matrix into a bml matrix.
 *
 * \ingroup convert_group
 *
 * \param N The number of rows/columns
 * \param matrix_precision The real precision
 * \param A The dense matrix
 * \return The bml matrix
 */
bml_matrix_ellblock_t
    * TYPED_FUNC(bml_import_from_dense_ellblock) (bml_dense_order_t order,
                                                  int N, void *A,
                                                  double threshold, int M,
                                                  bml_distribution_mode_t
                                                  distrib_mode)
{
    bml_matrix_ellblock_t *A_bml =
        TYPED_FUNC(bml_zero_matrix_ellblock) (N, M, distrib_mode);

    int NB = A_bml->NB;
    int MB = A_bml->MB;
    int *A_indexb = A_bml->indexb;
    int *A_nnzb = A_bml->nnzb;
    int *A_bsize = A_bml->bsize;

    REAL_T *dense_A = (REAL_T *) A;

    int *offset = malloc(NB * sizeof(int));
    offset[0] = 0;
    for (int ib = 1; ib < NB; ib++)
    {
        offset[ib] = offset[ib - 1] + A_bsize[ib - 1];
    }

    //allocate and assign data values
    for (int ib = 0; ib < NB; ib++)
    {
        A_nnzb[ib] = 0;
        for (int jb = 0; jb < NB; jb++)
        {
            int nelements = A_bsize[ib] * A_bsize[jb];
            REAL_T *A_ij = malloc(nelements * sizeof(REAL_T));
            switch (order)
            {
                case dense_row_major:
                    for (int ii = 0; ii < A_bsize[ib]; ii++)
                        for (int jj = 0; jj < A_bsize[jb]; jj++)
                            A_ij[ROWMAJOR(ii, jj, A_bsize[ib], A_bsize[jb])]
                                = dense_A[ROWMAJOR(offset[ib] + ii,
                                                   offset[jb] + jj, N, N)];
                    break;
                case dense_column_major:
                    for (int ii = 0; ii < A_bsize[ib]; ii++)
                        for (int jj = 0; jj < A_bsize[jb]; jj++)
                            A_ij[ROWMAJOR(ii, jj, A_bsize[ib], A_bsize[jb])]
                                = dense_A[COLMAJOR(offset[ib] + ii,
                                                   offset[jb] + jj, N, N)];
                    break;
                default:
                    LOG_ERROR("unknown order\n");
                    break;
            }
            double norminf = TYPED_FUNC(bml_norm_inf)
                (A_ij, A_bsize[ib], A_bsize[jb], A_bsize[jb]);
            if (norminf > threshold)
            {
                int ind = ROWMAJOR(ib, A_nnzb[ib], NB, MB);
                A_bml->ptr_value[ind] =
                    TYPED_FUNC(bml_allocate_block_ellblock) (A_bml, ib,
                                                             nelements);
                REAL_T *A_value = A_bml->ptr_value[ind];
                assert(A_value != NULL);
                memcpy(A_value, A_ij, nelements * sizeof(REAL_T));
                A_indexb[ind] = jb;
                A_nnzb[ib]++;
            }
            free(A_ij);
        }
    }
    free(offset);

    return A_bml;
}
