#include "../../macros.h"
#include "../../typed.h"
#include "../bml_allocate.h"
#include "../bml_normalize.h"
#include "../bml_parallel.h"
#include "../bml_types.h"
#include "bml_add_csr.h"
#include "bml_allocate_csr.h"
#include "bml_normalize_csr.h"
#include "bml_scale_csr.h"
#include "bml_types_csr.h"

#include <complex.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* Normalize csr matrix given Gershgorin bounds.
 *
 *  \ingroup normalize_group
 *
 *  \param A The matrix
 *  \param mineval Calculated min value
 *  \param maxeval Calculated max value
 */
void TYPED_FUNC(
    bml_normalize_csr) (
    bml_matrix_csr_t * A,
    double mineval,
    double maxeval)
{
    double maxminusmin = maxeval - mineval;
    double gershfact = maxeval / maxminusmin;
    REAL_T scalar = (REAL_T) - 1.0 / maxminusmin;
    double threshold = 0.0;

    bml_scale_inplace_csr(&scalar, A);
    bml_add_identity_csr(A, gershfact, threshold);
}

/** Calculate Gershgorin bounds for an csr matrix.
 *
 *  \ingroup normalize_group
 *
 *  \param A The matrix
 *  \param nrows Number of rows to use
 *  returns mineval Calculated min value
 *  returns maxeval Calculated max value
 */
void *TYPED_FUNC(
    bml_gershgorin_csr) (
    bml_matrix_csr_t * A)
{
    REAL_T radius, absham, dvalue;

    double emin = DBL_MAX;
    double emax = DBL_MIN;

    double *eval = bml_allocate_memory(sizeof(double) * 2);

    int N = A->N_;

    REAL_T rad[N];
    REAL_T dval[N];

#pragma omp parallel for                        \
  shared(N)         \
  shared(rad, dval)                             \
  private(absham, radius, dvalue)               \
  reduction(max:emax)                           \
  reduction(min:emin)
    for (int i = 0; i < N; i++)
    {
        radius = 0.0;
        dvalue = 0.0;
        int *cols = A->data_[i]->cols_;
        REAL_T *vals = (REAL_T *) A->data_[i]->vals_;
        const int annz = A->data_[i]->NNZ_;
        for (int pos = 0; pos < annz; pos++)
        {
            const int j = cols[pos];
            if (i == j)
                dvalue = vals[pos];
            else
            {
                absham = ABS(vals[pos]);
                radius += (double) absham;
            }
        }

        dval[i] = dvalue;
        rad[i] = radius;

/*
        emax =
            (emax >
             REAL_PART(dvalue + radius) ? emax : REAL_PART(dvalue + radius));
        emin =
            (emin <
             REAL_PART(dvalue - radius) ? emin : REAL_PART(dvalue - radius));

*/
    }
    for (int i = 0; i < N; i++)
    {
        if (REAL_PART(dval[i] + rad[i]) > emax)
            emax = REAL_PART(dval[i] + rad[i]);
        if (REAL_PART(dval[i] - rad[i]) < emin)
            emin = REAL_PART(dval[i] - rad[i]);
    }

    //printf("%d: emin = %e emax = %e\n", myRank, emin, emax);

    eval[0] = emin;
    eval[1] = emax;

    //printf("Global %d: emin = %e emax = %e\n", myRank, emin, emax);

    return eval;
}

/** Calculate Gershgorin bounds for a partial csr matrix.
 *
 *  \ingroup normalize_group
 *
 *  \param A The matrix
 *  \param nrows Number of rows to use
 *  returns mineval Calculated min value
 *  returns maxeval Calculated max value
 */
void *TYPED_FUNC(
    bml_gershgorin_partial_csr) (
    bml_matrix_csr_t * A,
    int nrows)
{
    REAL_T radius, absham, dvalue;

    double emin = DBL_MAX;
    double emax = DBL_MIN;

    double *eval = bml_allocate_memory(sizeof(double) * 2);

    int N = A->N_;

    REAL_T rad[N];
    REAL_T dval[N];

#pragma omp parallel for                        \
  shared(N)         \
  shared(rad, dval)                             \
  private(absham, radius, dvalue)               \
  reduction(max:emax)                           \
  reduction(min:emin)
    for (int i = 0; i < nrows; i++)
    {
        radius = 0.0;
        dvalue = 0.0;
        int *cols = A->data_[i]->cols_;
        REAL_T *vals = (REAL_T *) A->data_[i]->vals_;
        const int annz = A->data_[i]->NNZ_;
        for (int pos = 0; pos < annz; pos++)
        {
            const int j = cols[pos];
            if (i == j)
                dvalue = vals[pos];
            else
            {
                absham = ABS(vals[pos]);
                radius += (double) absham;
            }

        }

        dval[i] = dvalue;
        rad[i] = radius;

    }

    for (int i = 0; i < nrows; i++)
    {
        if (REAL_PART(dval[i] + rad[i]) > emax)
            emax = REAL_PART(dval[i] + rad[i]);
        if (REAL_PART(dval[i] - rad[i]) < emin)
            emin = REAL_PART(dval[i] - rad[i]);

    }

    eval[0] = emin;
    eval[1] = emax;

    return eval;
}
