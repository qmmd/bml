#include "../../macros.h"
#include "../../typed.h"
#include "../blas.h"
#include "../bml_allocate.h"
#include "../bml_logger.h"
#include "../bml_parallel.h"
#include "../bml_scale.h"
#include "../bml_types.h"
#include "bml_allocate_ellblock.h"
#include "bml_copy_ellblock.h"
#include "bml_scale_ellblock.h"
#include "bml_types_ellblock.h"

#include <stdlib.h>
#include <string.h>
#include <assert.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/** Scale an ellblock matrix - result is a new matrix.
 *
 *  \ingroup scale_group
 *
 *  \param A The matrix to be scaled
 *  \return A scale version of matrix A.
 */
bml_matrix_ellblock_t *TYPED_FUNC(
    bml_scale_ellblock_new) (
    void *_scale_factor,
    bml_matrix_ellblock_t * A)
{
    REAL_T *scale_factor = _scale_factor;
    bml_matrix_ellblock_t *B = TYPED_FUNC(bml_copy_ellblock_new) (A);

    TYPED_FUNC(bml_scale_ellblock) (scale_factor, A, B);

    return B;
}

/** Scale an ellblock matrix.
 *
 *  \ingroup scale_group
 *
 *  \param A The matrix to be scaled
 *  \param B Scaled version of matrix A
 */
void TYPED_FUNC(
    bml_scale_ellblock) (
    void *_scale_factor,
    bml_matrix_ellblock_t * A,
    bml_matrix_ellblock_t * B)
{
#ifdef NOBLAS
    LOG_ERROR("No BLAS library");
#else
    if (A != B)
    {
        TYPED_FUNC(bml_copy_ellblock) (A, B);
    }

    REAL_T *scale_factor = _scale_factor;
    REAL_T **B_ptr_value = (REAL_T **) B->ptr_value;
    int inc = 1;

    for (int ib = 0; ib < A->NB; ib++)
    {
        for (int jp = 0; jp < A->nnzb[ib]; jp++)
        {
            int ind = ROWMAJOR(ib, jp, A->NB, A->MB);
            int jb = A->indexb[ind];
            int nElems = A->bsize[ib] * A->bsize[jb];

            REAL_T *B_value = B_ptr_value[ind];
            assert(B_value != NULL);

            C_BLAS(SCAL) (&nElems, scale_factor, &(B_value[0]), &inc);
        }
    }
#endif
}

void TYPED_FUNC(
    bml_scale_inplace_ellblock) (
    void *scale_factor,
    bml_matrix_ellblock_t * A)
{
    TYPED_FUNC(bml_scale_ellblock) (scale_factor, A, A);
}
