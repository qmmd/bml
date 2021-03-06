#include "../../macros.h"
#include "../../typed.h"
#include "../bml_allocate.h"
#include "../bml_parallel.h"
#include "../bml_transpose.h"
#include "../bml_types.h"
#include "bml_allocate_ellpack.h"
#include "bml_transpose_ellpack.h"
#include "bml_types_ellpack.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define COMPUTE_ON_HOST

/** Transpose a matrix.
 *
 *  \ingroup transpose_group
 *
 *  \param A The matrix to be transposed
 *  \return the transposed A
 */
bml_matrix_ellpack_t
    * TYPED_FUNC(bml_transpose_new_ellpack) (bml_matrix_ellpack_t * A)
{
    bml_matrix_dimension_t matrix_dimension = { A->N, A->N, A->M };

    bml_matrix_ellpack_t *B =
        TYPED_FUNC(bml_noinit_matrix_ellpack) (matrix_dimension,
                                               A->distribution_mode);

    REAL_T *A_value = (REAL_T *) A->value;
    int *A_index = A->index;
    int *A_nnz = A->nnz;
    int *A_localRowMin = A->domain->localRowMin;
    int *A_localRowMax = A->domain->localRowMax;

    REAL_T *B_value = (REAL_T *) B->value;
    int *B_index = B->index;
    int *B_nnz = B->nnz;

    int N = A->N;
    int M = A->M;

    int myRank = bml_getMyRank();

    // Transpose all elements
#ifdef _OPENMP
    omp_lock_t *row_lock = (omp_lock_t *) malloc(sizeof(omp_lock_t) * N);

#pragma omp parallel for
    for (int i = 0; i < N; i++)
    {
        omp_init_lock(&row_lock[i]);
    }
#endif

#ifdef USE_OMP_OFFLOAD
#ifdef COMPUTE_ON_HOST
#pragma omp target update from(A_index[:N*M], A_value[:N*M], A_nnz[:N])
#else
#pragma omp target map(to:row_lock[:N])
#endif
#endif
    {                           // begin target region
#pragma omp parallel for                                \
   shared(B_index, B_value, B_nnz)                      \
   shared(A_index, A_value, A_nnz, row_lock)
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < A_nnz[i]; j++)
            {
                int trow = A_index[ROWMAJOR(i, j, N, M)];
#ifdef _OPENMP
                omp_set_lock(&row_lock[trow]);
#endif
                int colcnt = B_nnz[trow];
                B_index[ROWMAJOR(trow, colcnt, N, M)] = i;
                B_value[ROWMAJOR(trow, colcnt, N, M)] =
                    A_value[ROWMAJOR(i, j, N, M)];
                B_nnz[trow]++;
#ifdef _OPENMP
                omp_unset_lock(&row_lock[trow]);
#endif
            }
        }
    }                           // end target region

#if defined (USE_OMP_OFFLOAD) && defined(COMPUTE_ON_HOST)
#pragma omp target update to(B_index[:N*M], B_value[:N*M], B_nnz[:N])
#endif
    return B;
    /*
       int Alrmin = A_localRowMin[myRank];
       int Alrmax = A_localRowMax[myRank];

       #pragma omp parallel for               \
       shared(N, M, B_index, B_value, B_nnz) \
       shared(A_index, A_value, A_nnz,Alrmin,Alrmax)
       //for (int i = 0; i < N; i++)

       for (int i = Alrmin; i < Alrmax; i++)
       {
       for (int j = 0; j < N; j++)
       {
       int Annzj = A_nnz[j];
       for (int k = 0; k < Annzj; k++)
       {
       if (A_index[ROWMAJOR(j, k, N, M)] != i) {}
       else {
       B_index[ROWMAJOR(i, B_nnz[i], N, M)] = j;
       B_value[ROWMAJOR(i, B_nnz[i], N, M)] = A_value[ROWMAJOR(j, k, N, M)];
       B_nnz[i]++;
       break;
       }
       }
       }
       }

       return B;
     */
}


/** Transpose a matrix in place.
 *
 *  \ingroup transpose_group
 *
 *  \param A The matrix to be transposeed
 *  \return the transposed A
 */
void TYPED_FUNC(
    bml_transpose_ellpack) (
    bml_matrix_ellpack_t * A)
{
    int N = A->N;
    int M = A->M;

    REAL_T *A_value = (REAL_T *) A->value;
    int *A_index = A->index;
    int *A_nnz = A->nnz;

#if defined(USE_OMP_OFFLOAD)
#ifdef COMPUTE_ON_HOST
#pragma omp target update from(A_index[:N*M], A_value[:N*M], A_nnz[:N])
#else
#pragma omp target
#endif
#endif
    {                           // begin target region
#pragma omp parallel for shared(N, M, A_value, A_index, A_nnz)
        for (int i = 0; i < N; i++)
        {
            for (int j = A_nnz[i] - 1; j >= 0; j--)
            {
                if (A_index[ROWMAJOR(i, j, N, M)] > i)
                {
                    int ind = A_index[ROWMAJOR(i, j, N, M)];
                    int exchangeDone = 0;
                    for (int k = 0; k < A_nnz[ind]; k++)
                    {
                        // Existing corresponding value for transpose - exchange
                        if (A_index[ROWMAJOR(ind, k, N, M)] == i)
                        {
                            REAL_T tmp = A_value[ROWMAJOR(i, j, N, M)];

#pragma omp critical
                            {
                                A_value[ROWMAJOR(i, j, N, M)] =
                                    A_value[ROWMAJOR(ind, k, N, M)];
                                A_value[ROWMAJOR(ind, k, N, M)] = tmp;
                            }
                            exchangeDone = 1;
                            break;
                        }
                    }

                    // If no match add to end of row
                    if (!exchangeDone)
                    {
                        int jind = A_nnz[ind];

#pragma omp critical
                        {
                            A_index[ROWMAJOR(ind, jind, N, M)] = i;
                            A_value[ROWMAJOR(ind, jind, N, M)] =
                                A_value[ROWMAJOR(i, j, N, M)];
                            A_nnz[ind]++;
                            A_nnz[i]--;
                        }
                    }
                }
            }
        }
    }                           // end target region
#ifdef USE_OMP_OFFLOAD
#ifdef COMPUTE_ON_HOST
#pragma omp target update to(A_index[:N*M], A_value[:N*M], A_nnz[:N])
#endif
#endif

}
