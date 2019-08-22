#include "../../macros.h"
#include "../../typed.h"
#include "../bml_allocate.h"
#include "../bml_logger.h"
#include "../bml_types.h"
#include "bml_allocate_ellsort.h"
#include "bml_import_ellsort.h"
#include "bml_types_ellsort.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/** Convert a dense matrix into a bml matrix.
 *
 * \ingroup convert_group
 *
 * \param N The number of rows/columns
 * \param matrix_precision The real precision
 * \param A The dense matrix
 * \return The bml matrix
 */
bml_matrix_ellsort_t *TYPED_FUNC(
    bml_import_from_dense_ellsort) (
    bml_dense_order_t order,
    int N,
    void *A,
    double threshold,
    int M,
    bml_distribution_mode_t distrib_mode)
{
    bml_matrix_ellsort_t *A_bml =
        TYPED_FUNC(bml_zero_matrix_ellsort) (N, M, distrib_mode);

    int *A_index = A_bml->index;
    int *A_nnz = A_bml->nnz;

    REAL_T *dense_A = (REAL_T *) A;
    REAL_T *A_value = A_bml->value;

#pragma omp parallel for shared(A_value, A_index, A_nnz, dense_A)
    for (int i = 0; i < N; i++)
    {
        A_nnz[i] = 0;
        for (int j = 0; j < N; j++)
        {
            REAL_T A_ij;
            switch (order)
            {
                case dense_row_major:
                    A_ij = dense_A[ROWMAJOR(i, j, N, N)];
                    break;
                case dense_column_major:
                    A_ij = dense_A[COLMAJOR(i, j, N, N)];
                    break;
                default:
                    LOG_ERROR("unknown order\n");
                    break;
            }
            if (is_above_threshold(A_ij, threshold))
            {
                A_value[ROWMAJOR(i, A_nnz[i], N, M)] = A_ij;
                A_index[ROWMAJOR(i, A_nnz[i], N, M)] = j;
                A_nnz[i]++;
            }
        }
    }

    return A_bml;
}
