#include "../../macros.h"
#include "../../typed.h"
#include "../bml_introspection.h"
#include "../bml_logger.h"
#include "bml_adjungate_triangle_ellpack.h"
#include "bml_types.h"
#include "bml_types_ellpack.h"

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>


#ifdef _OPENMP
#include <omp.h>
#endif

/** Adjungates a triangle of a matrix in place.
 *
 *  \ingroup adjungate_triangle_group
 *
 *  \param A  The matrix for which the triangle should be adjungated
 *  \param triangle  Which triangle to adjungate ('u': upper, 'l': lower)
 *  WARNING: Please verify race conditions and parallel performances ....
 */
void TYPED_FUNC(
    bml_adjungate_triangle_ellpack) (
    bml_matrix_ellpack_t * A,
    char *triangle)
{
    int A_N = A->N;
    int A_M = A->M;
    int l;
    int ll;
    int j;
    REAL_T *A_value = (REAL_T *) A->value;
    int *A_index = A->index;
    int *A_nnz = A->nnz;
#ifdef _OPENMP
    omp_lock_t lock[A_M];
#endif

    switch (*triangle)
    {
        case 'u':
#ifdef _OPENMP
            for (int i = 0; i < A_M; i++)
                omp_init_lock(&(lock[i]));
#endif

#pragma omp parallel for default(none) shared(A_N,A_M,A_index,A_nnz,A_value,lock) \
      private(j,l,ll)
            //    WARNING: Please, check for race conditions ...

            for (int i = 0; i < A_N; i++)       // For every row
            {
                l = A_nnz[i];
                for (int j = 0; j < l; j++)     // We search for indices gt 0.
                {
                    ll = A_index[ROWMAJOR(i, j, A_N, A_M)];
                    if (ll > 0)
                    {
                        if (ll > i)
                        {
#ifdef _OPENMP
                            omp_set_lock(&(lock[ll]));
#endif
                            A_index[ROWMAJOR(ll, A_nnz[ll], A_N, A_M)] = i;
                            A_value[ROWMAJOR(ll, A_nnz[ll], A_N, A_M)] =
                                conj(A_value[ROWMAJOR(i, j, A_N, A_M)]);
                            A_nnz[ll]++;
#ifdef _OPENMP
                            omp_unset_lock(&(lock[ll]));
#endif

                        }
                    }
                }
            }
#ifdef _OPENMP
            for (int i = 0; i < A_M; i++)
                omp_destroy_lock(&(lock[i]));
#endif

            break;

        case 'l':
#ifdef _OPENMP
            for (int i = 0; i < A_M; i++)
                omp_init_lock(&(lock[i]));
#endif

#pragma omp parallel for default(none) shared(lock,A_N,A_M,A_index,A_nnz,A_value) \
      private(j,l,ll)
            //    WARNING: Please, check for race conditions and parallel performances ...
            for (int i = 0; i < A_N; i++)
            {
                l = A_nnz[i];
                for (int j = 0; j < l; j++)
                {
                    ll = A_index[ROWMAJOR(i, j, A_N, A_M)];
                    if (ll >= 0)
                    {
                        if (ll < i)
                        {
#ifdef _OPENMP
                            omp_set_lock(&(lock[ll]));
#endif
                            A_index[ROWMAJOR(ll, A_nnz[ll], A_N, A_M)] = i;
                            A_value[ROWMAJOR(ll, A_nnz[ll], A_N, A_M)] =
                                conj(A_value[ROWMAJOR(i, j, A_N, A_M)]);
                            A_nnz[ll]++;
#ifdef _OPENMP
                            omp_unset_lock(&(lock[ll]));
#endif

                        }
                    }
                }
            }
#ifdef _OPENMP
            for (int i = 0; i < A_M; i++)
                omp_destroy_lock(&(lock[i]));
#endif

            break;
        default:
            LOG_ERROR("unknown triangle %c\n", triangle);
            break;
    }
}
