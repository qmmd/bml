#include "../../macros.h"
#include "../../typed.h"
#include "../bml_logger.h"
#include "../bml_parallel.h"
#include "../bml_submatrix.h"
#include "../bml_trace.h"
#include "../bml_types.h"
#include "bml_trace_ellblock.h"
#include "bml_types_ellblock.h"

#include <complex.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/** Calculate the trace of a matrix.
 *
 *  \ingroup trace_group
 *
 *  \param A The matrix to calculate a trace for
 *  \return the trace of A
 */
double TYPED_FUNC(
    bml_trace_ellblock) (
    bml_matrix_ellblock_t * A)
{
    int NB = A->NB;
    int MB = A->MB;

    int *A_indexb = (int *) A->indexb;
    int *A_nnzb = (int *) A->nnzb;
    int *bsize = A->bsize;

    REAL_T trace = 0.0;
    REAL_T **A_ptr_value = (REAL_T **) A->ptr_value;

    for (int ib = 0; ib < NB; ib++)
    {
        for (int jp = 0; jp < A_nnzb[ib]; jp++)
        {
            int ind = ROWMAJOR(ib, jp, NB, MB);
            int jb = A_indexb[ind];
            if (ib == jb)
            {
                REAL_T *A_value = A_ptr_value[ind];
                for (int ii = 0; ii < bsize[ib]; ii++)
                    trace += A_value[ROWMAJOR(ii, ii, bsize[ib], bsize[ib])];
                break;
            }
        }
    }

    return (double) REAL_PART(trace);
}

/** Calculate the trace of a matrix multiplication.
 * Both matrices must have the same size.
 *
 *  \ingroup trace_group
 *
 *  \param A The matrix A
 *  \param A The matrix B
 *  \return the trace of A*B
 */
double TYPED_FUNC(
    bml_trace_mult_ellblock) (
    bml_matrix_ellblock_t * A,
    bml_matrix_ellblock_t * B)
{
    int NB = A->NB;
    int MB = A->MB;

    int *A_indexb = (int *) A->indexb;
    int *A_nnzb = (int *) A->nnzb;
    int *bsize = A->bsize;

    int *B_indexb = (int *) B->indexb;
    int *B_nnzb = (int *) B->nnzb;

    REAL_T trace = 0.0;
    REAL_T **A_ptr_value = (REAL_T **) A->ptr_value;
    REAL_T **B_ptr_value = (REAL_T **) B->ptr_value;

    if (NB != B->NB || MB != B->MB)
    {
        LOG_ERROR
            ("bml_trace_mult_ellblock: Matrices A and B have different number of blocks.");
    }

    for (int ib = 0; ib < NB; ib++)
    {
        for (int jp = 0; jp < A_nnzb[ib]; jp++)
        {
            int indA = ROWMAJOR(ib, jp, NB, MB);
            int jb = A_indexb[indA];
            REAL_T *A_value = A_ptr_value[indA];

            for (int kp = 0; kp < B_nnzb[jb]; kp++)
            {
                int indB = ROWMAJOR(jb, kp, NB, MB);
                int kb = B_indexb[indB];
                if (kb == ib)
                {
                    REAL_T *B_value = B_ptr_value[indB];
                    for (int ii = 0; ii < bsize[ib]; ii++)
                        for (int jj = 0; jj < bsize[jb]; jj++)
                        {
                            int ka = ROWMAJOR(ii, jj, bsize[ib], bsize[jb]);
                            int kb = ROWMAJOR(jj, ii, bsize[jb], bsize[ib]);
                            trace += A_value[ka] * B_value[kb];
                        }
                }
            }
        }
    }

    return (double) REAL_PART(trace);
}
