#include "../../macros.h"
#include "../../typed.h"
#include "../bml_introspection.h"
#include "../bml_types.h"
#include "bml_introspection_ellblock.h"
#include "bml_types_ellblock.h"
#include "bml_utilities_ellblock.h"

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

/** Return the sparsity of a matrix.
 *
 *  Note that the the sparsity of a matrix is defined
 *  as NumberOfZeroes/N*N where N is the matrix dimension.
 *  To density of matrix A will be defined as 1-sparsity(A)
 *
 * \ingroup introspection_group_C
 *
 * \param A The bml matrix.
 * \param threshold The threshold used to compute the sparsity.
 * \return The sparsity of A.
 */
double TYPED_FUNC(
    bml_get_sparsity_ellblock) (
    bml_matrix_ellblock_t * A,
    double threshold)
{
    int nnzs = 0;
    double sparsity;
    REAL_T **A_ptr_value = (REAL_T **) A->ptr_value;
    int NB = A->NB;
    int MB = A->MB;
    int *A_nnzb = A->nnzb;
    int *A_indexb = A->indexb;
    int *bsize = A->bsize;

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
                    int index = ROWMAJOR(ii, jj, bsize[ib], bsize[jb]);
                    if (ABS(A_value[index]) > threshold)
                    {
                        nnzs++;
                    }
                }
        }
    }

    sparsity = (1.0 - (double) nnzs / ((double) (A->N * A->N)));
    return sparsity;
}
