#include "../../macros.h"
#include "../../typed.h"
#include "../bml_add.h"
#include "../bml_allocate.h"
#include "../bml_normalize.h"
#include "../bml_parallel.h"
#include "../bml_types.h"
#include "bml_add_dense.h"
#include "bml_allocate_dense.h"
#include "bml_getters_dense.h"
#include "bml_normalize_dense.h"
#include "bml_scale_dense.h"
#include "bml_types_dense.h"

#ifdef BML_USE_MAGMA
#include "magma_v2.h"
#endif

#include <complex.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/** Normalize dense matrix given Gershgorin bounds.
 *
 *  \ingroup normalize_group
 *
 *  \param A The matrix
 *  \param mineval Calculated min value
 *  \param maxeval Calculated max value
 */
void TYPED_FUNC(
    bml_normalize_dense) (
    bml_matrix_dense_t * A,
    double mineval,
    double maxeval)
{
    double maxminusmin = maxeval - mineval;
    double beta = maxeval / maxminusmin;
    double alpha = (double) -1.0 / maxminusmin;

    TYPED_FUNC(bml_scale_add_identity_dense) (A, alpha, beta);
}

/** Calculate Gershgorin bounds for a dense matrix.
 *
 *  \ingroup gershgorin_group
 *
 *  \param A The matrix
 *  returns mineval Calculated min value
 *  returns maxeval Calculated max value
 */
void *TYPED_FUNC(
    bml_gershgorin_dense) (
    bml_matrix_dense_t * A)
{
    REAL_T radius, dvalue, absham;

    int N = A->N;

    int *A_localRowMin = A->domain->localRowMin;
    int *A_localRowMax = A->domain->localRowMax;

    int myRank = bml_getMyRank();

    double emin = DBL_MAX;
    double emax = DBL_MIN;

    double *eval = bml_allocate_memory(sizeof(double) * 2);
#ifdef BML_USE_MAGMA
    //copy data from GPU to CPU to do the work on the CPU
    REAL_T *A_matrix = bml_allocate_memory(sizeof(REAL_T) * A->N * A->N);

    MAGMA(getmatrix) (A->N, A->N,
                      A->matrix, A->ld, (MAGMA_T *) A_matrix, A->N, A->queue);
#else
    REAL_T *A_matrix = A->matrix;
#endif

#pragma omp parallel for                        \
  shared(N, A_matrix)                           \
  shared(A_localRowMin, A_localRowMax, myRank)  \
  private(absham, radius, dvalue)               \
  reduction(max:emax)                           \
  reduction(min:emin)
    //for (int i = 0; i < N; i++)
    for (int i = A_localRowMin[myRank]; i < A_localRowMax[myRank]; i++)
    {
        radius = 0.0;

        for (int j = 0; j < N; j++)
        {
            absham = ABS(A_matrix[ROWMAJOR(i, j, N, N)]);
            radius += (double) absham;
        }

        dvalue = A_matrix[ROWMAJOR(i, i, N, N)];

        radius -= ABS(dvalue);

/*
        emax =
            (emax >
             REAL_PART(dvalue + radius) ? emax : REAL_PART(dvalue + radius));
        emin =
            (emin <
             REAL_PART(dvalue - radius) ? emin : REAL_PART(dvalue - radius));
*/
        if (REAL_PART(dvalue + radius) > emax)
            emax = REAL_PART(dvalue + radius);
        if (REAL_PART(dvalue - radius) < emin)
            emin = REAL_PART(dvalue - radius);
    }

#ifdef BML_USE_MAGMA
    bml_free_memory(A_matrix);
#endif

#ifdef DO_MPI
    if (bml_getNRanks() > 1 && A->distribution_mode == distributed)
    {
        bml_minRealReduce(&emin);
        bml_maxRealReduce(&emax);
    }
#endif

    eval[0] = emin;
    eval[1] = emax;

    return eval;
}

/** Calculate Gershgorin bounds for a partial dense matrix.
 *
 *  \ingroup gershgorin_group
 *
 *  \param A The matrix
 *  \param nrows Number of rows used
 *  returns mineval Calculated min value
 *  returns maxeval Calculated max value
 */
void *TYPED_FUNC(
    bml_gershgorin_partial_dense) (
    bml_matrix_dense_t * A,
    int nrows)
{
    REAL_T radius, dvalue, absham;

    int N = A->N;
    REAL_T *A_matrix = A->matrix;

    double emin = DBL_MAX;
    double emax = DBL_MIN;

    double *eval = bml_allocate_memory(sizeof(double) * 2);

#pragma omp parallel for                        \
  shared(N, A_matrix)                           \
  private(absham, radius, dvalue)               \
  reduction(max:emax)                           \
  reduction(min:emin)
    for (int i = 0; i < nrows; i++)
    {
        radius = 0.0;

        for (int j = 0; j < N; j++)
        {
            absham = ABS(A_matrix[ROWMAJOR(i, j, N, N)]);
            radius += (double) absham;
        }

        dvalue = A_matrix[ROWMAJOR(i, i, N, N)];

        radius -= ABS(dvalue);

        if (REAL_PART(dvalue + radius) > emax)
            emax = REAL_PART(dvalue + radius);
        if (REAL_PART(dvalue - radius) < emin)
            emin = REAL_PART(dvalue - radius);
    }

    eval[0] = emin;
    eval[1] = emax;

    return eval;
}
