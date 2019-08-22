#include "../../macros.h"
#include "../../typed.h"
#include "../bml_adjungate_triangle.h"
#include "../bml_allocate.h"
#include "../bml_logger.h"
#include "../bml_types.h"
#include "bml_adjungate_triangle_dense.h"
#include "bml_allocate_dense.h"
#include "bml_types_dense.h"

#include <complex.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/** Adjungates a triangle of a matrix in place.
 *
 *  \ingroup adjungate_triangle_group
 *
 *  \param A[in,out]  The matrix for which the triangle should be adjungated
 *  \param triangle[in]  Which triangle to adjungate ('u': upper, 'l': lower)
 */
void TYPED_FUNC(
    bml_adjungate_triangle_dense) (
    bml_matrix_dense_t * A,
    char *triangle)
{
    int N = A->N;

    REAL_T *A_matrix = A->matrix;

    switch (*triangle)
    {
        case 'u':
#pragma omp parallel for shared(N, A_matrix)
            for (int i = 0; i < N - 1; i++)
            {
                for (int j = i + 1; j < N; j++)
                {
                    A_matrix[ROWMAJOR(j, i, N, N)] =
                        conj(A_matrix[ROWMAJOR(i, j, N, N)]);
                }
            }
            break;

        case 'l':
#pragma omp parallel for shared(N, A_matrix)
            for (int i = 0; i < N - 1; i++)
            {
                for (int j = i + 1; j < N; j++)
                {
                    A_matrix[ROWMAJOR(i, j, N, N)] =
                        conj(A_matrix[ROWMAJOR(j, i, N, N)]);
                }
            }
            break;

        default:
            LOG_ERROR("Unknown triangle type in bml_adjungate\n");
    }
}
