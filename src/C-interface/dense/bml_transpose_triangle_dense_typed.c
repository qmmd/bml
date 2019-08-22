#include "../../macros.h"
#include "../../typed.h"
#include "../bml_allocate.h"
#include "../bml_transpose_triangle.h"
#include "../bml_types.h"
#include "bml_allocate_dense.h"
#include "bml_transpose_triangle_dense.h"
#include "bml_types_dense.h"

#include <complex.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/** Transposes a triangle of a matrix in place.
 *
 *  \ingroup transpose_triangle_group
 *
 *  \param A  The matrix for which the triangle should be transposed
 *  \param triangle  Which triangle to transpose ('u': upper, 'l': lower)
 */
void TYPED_FUNC(
    bml_transpose_triangle_dense) (
    bml_matrix_dense_t * A,
    char triangle)
{
    int N = A->N;

    REAL_T *A_matrix = A->matrix;

    switch (triangle)
    {
        case 'u':
#pragma omp parallel for shared(N, A_matrix)
            for (int i = 0; i < N - 1; i++)
            {
                for (int j = i + 1; j < N; j++)
                {
                    A_matrix[ROWMAJOR(j, i, N, N)] =
                        A_matrix[ROWMAJOR(i, j, N, N)];
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
                        A_matrix[ROWMAJOR(j, i, N, N)];
                }
            }
            break;
    }
}
