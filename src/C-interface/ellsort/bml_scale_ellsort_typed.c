#include "../../typed.h"
#include "../blas.h"
#include "../bml_allocate.h"
#include "../bml_logger.h"
#include "../bml_parallel.h"
#include "../bml_scale.h"
#include "../bml_types.h"
#include "bml_allocate_ellsort.h"
#include "bml_copy_ellsort.h"
#include "bml_scale_ellsort.h"
#include "bml_types_ellsort.h"

#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/** Scale an ellsort matrix - result is a new matrix.
 *
 *  \ingroup scale_group
 *
 *  \param A The matrix to be scaled
 *  \return A scale version of matrix A.
 */
bml_matrix_ellsort_t *TYPED_FUNC(
    bml_scale_ellsort_new) (
    void *_scale_factor,
    bml_matrix_ellsort_t * A)
{
    REAL_T *scale_factor = _scale_factor;
    bml_matrix_ellsort_t *B = TYPED_FUNC(bml_copy_ellsort_new) (A);

    REAL_T *B_value = B->value;
    int myRank = bml_getMyRank();
    int nElems = B->domain->localRowExtent[myRank] * B->M;
    int startIndex = B->domain->localDispl[myRank];
    int inc = 1;

#ifdef NOBLAS
    LOG_ERROR("No BLAS library");
#else
    C_BLAS(SCAL) (&nElems, scale_factor, &(B_value[startIndex]), &inc);
#endif

    return B;
}

/** Scale an ellsort matrix.
 *
 *  \ingroup scale_group
 *
 *  \param A The matrix to be scaled
 *  \param B Scaled version of matrix A
 */
void TYPED_FUNC(
    bml_scale_ellsort) (
    void *_scale_factor,
    bml_matrix_ellsort_t * A,
    bml_matrix_ellsort_t * B)
{
    if (A != B)
    {
        TYPED_FUNC(bml_copy_ellsort) (A, B);
    }

    REAL_T *scale_factor = _scale_factor;
    REAL_T *B_value = B->value;
    int myRank = bml_getMyRank();
    int nElems = B->domain->localRowExtent[myRank] * B->M;
    int startIndex = B->domain->localDispl[myRank];
    int inc = 1;

#ifdef NOBLAS
    LOG_ERROR("No BLAS library");
#else
    C_BLAS(SCAL) (&nElems, scale_factor, &(B_value[startIndex]), &inc);
#endif

}

void TYPED_FUNC(
    bml_scale_inplace_ellsort) (
    void *_scale_factor,
    bml_matrix_ellsort_t * A)
{
    REAL_T *scale_factor = _scale_factor;
    REAL_T *A_value = A->value;

    int myRank = bml_getMyRank();
    int number_elements = A->domain->localRowExtent[myRank] * A->M;
    int startIndex = A->domain->localDispl[myRank];
    int inc = 1;

#ifdef NOBLAS
    LOG_ERROR("No BLAS library");
#else
    C_BLAS(SCAL) (&number_elements, scale_factor, &(A_value[startIndex]),
                  &inc);
#endif

}
