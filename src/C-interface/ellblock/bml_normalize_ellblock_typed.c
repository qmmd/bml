#include "../../macros.h"
#include "../../typed.h"
#include "../bml_logger.h"
#include "../bml_allocate.h"
#include "../bml_normalize.h"
#include "../bml_parallel.h"
#include "../bml_types.h"
#include "bml_allocate_ellblock.h"
#include "bml_normalize_ellblock.h"
#include "bml_scale_ellblock.h"
#include "bml_add_ellblock.h"
#include "bml_types_ellblock.h"

#include <complex.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* Normalize ellblock matrix given Gershgorin bounds.
 *
 *  \ingroup normalize_group
 *
 *  \param A The matrix
 *  \param mineval Calculated min value
 *  \param maxeval Calculated max value
 */
void TYPED_FUNC(
    bml_normalize_ellblock) (
    bml_matrix_ellblock_t * A,
    double mineval,
    double maxeval)
{
    double maxminusmin = maxeval - mineval;
    double gershfact = maxeval / maxminusmin;
    REAL_T scalar = (REAL_T) - 1.0 / maxminusmin;
    double threshold = 0.0;

    bml_scale_inplace_ellblock(&scalar, A);
    bml_add_identity_ellblock(A, gershfact, threshold);
}

/** Calculate Gershgorin bounds for an ellblock matrix.
 *
 *  \ingroup normalize_group
 *
 *  \param A The matrix
 *  \param nrows Number of rows to use
 *  returns mineval Calculated min value
 *  returns maxeval Calculated max value
 */
void *TYPED_FUNC(
    bml_gershgorin_ellblock) (
    bml_matrix_ellblock_t * A)
{
    double emin = DBL_MAX;
    double emax = DBL_MIN;

    double *eval = bml_allocate_memory(sizeof(double) * 2);

    int NB = A->NB;
    int MB = A->MB;
    int *A_nnzb = (int *) A->nnzb;
    int *A_indexb = (int *) A->indexb;
    int *bsize = A->bsize;

    REAL_T *rad = calloc(A->N, sizeof(REAL_T));
    REAL_T *dval = calloc(A->N, sizeof(REAL_T));

    REAL_T **A_ptr_value = (REAL_T **) A->ptr_value;

    int ioffset = 0;
    for (int ib = 0; ib < NB; ib++)
    {
        for (int jp = 0; jp < A_nnzb[ib]; jp++)
        {
            int ind = ROWMAJOR(ib, jp, NB, MB);
            int jb = A_indexb[ind];
            REAL_T *A_value = A_ptr_value[ind];
            for (int ii = 0; ii < bsize[ib]; ii++)
            {
                int i = ioffset + ii;
                for (int jj = 0; jj < bsize[jb]; jj++)
                {
                    REAL_T absham =
                        A_value[ROWMAJOR(ii, jj, bsize[ib], bsize[jb])];
                    if (ib == jb && ii == jj)
                        dval[i] = (double) absham;
                    else
                        rad[i] += (double) ABS(absham);
                }
            }
        }

        ioffset += bsize[ib];
    }
    for (int i = 0; i < A->N; i++)
    {
        if (REAL_PART(dval[i] + rad[i]) > emax)
            emax = REAL_PART(dval[i] + rad[i]);
        if (REAL_PART(dval[i] - rad[i]) < emin)
            emin = REAL_PART(dval[i] - rad[i]);

    }

    //printf("%d: emin = %e emax = %e\n", myRank, emin, emax);

#ifdef DO_MPI
    if (bml_getNRanks() > 1 && A->distribution_mode == distributed)
    {
        bml_minRealReduce(&emin);
        bml_maxRealReduce(&emax);
    }
#endif

    eval[0] = emin;
    eval[1] = emax;

    //printf("Global %d: emin = %e emax = %e\n", myRank, emin, emax);
    free(dval);
    free(rad);

    return eval;
}
