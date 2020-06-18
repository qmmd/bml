#include "../../macros.h"
#include "../../typed.h"
#include "../bml_allocate.h"
#include "../bml_getters.h"
#include "../bml_introspection.h"
#include "../bml_logger.h"
#include "../bml_setters.h"
#include "bml_allocate_ellblock.h"
#include "bml_setters_ellblock.h"
#include "bml_types_ellblock.h"

#include <complex.h>
#include <assert.h>
#include <stdio.h>

bml_matrix_ellblock_t *TYPED_FUNC(
    bml_convert_ellblock) (
    bml_matrix_t * A,
    bml_matrix_precision_t matrix_precision,
    int M,
    bml_distribution_mode_t distrib_mode)
{
    int N = bml_get_N(A);

    if (N < 0)
    {
        LOG_ERROR("A is not intialized\n");
    }

    //create new empty matrix
    int *bsize = bml_get_block_sizes(N, N);
    int nb = bml_get_nb();

    bml_matrix_ellblock_t *B =
        TYPED_FUNC(bml_block_matrix_ellblock) (nb, nb, M, bsize,
                                               distrib_mode);

    int NB = B->NB;
    int *offset = malloc(NB * sizeof(int));
    offset[0] = 0;
    for (int ib = 1; ib < NB; ib++)
        offset[ib] = offset[ib - 1] + bsize[ib - 1];

    REAL_T **B_ptr_value = (REAL_T **) B->ptr_value;

    for (int ib = 0; ib < NB; ib++)
    {
        B->nnzb[ib] = NB;
        for (int jb = 0; jb < NB; jb++)
        {
            int nelements = bsize[ib] * bsize[jb];
            int ind = ROWMAJOR(ib, jb, NB, NB);
            B_ptr_value[ind]
                = TYPED_FUNC(bml_allocate_block_ellblock) (B, ib, nelements);
            B->indexb[ind] = jb;

            REAL_T *block = malloc(nelements * sizeof(REAL_T));
            for (int ii = 0; ii < bsize[ib]; ii++)
                for (int jj = 0; jj < bsize[jb]; jj++)
                {
                    int i = offset[ib] + ii;
                    int j = offset[jb] + jj;
                    REAL_T alpha = *(REAL_T *) bml_get_element(A, i, j);
                    block[ROWMAJOR(ii, jj, bsize[ib], bsize[jb])] = alpha;
                }
            TYPED_FUNC(bml_set_block_ellblock) (B, ib, jb, block);
            free(block);
        }
    }
    free(offset);

    return B;
}
