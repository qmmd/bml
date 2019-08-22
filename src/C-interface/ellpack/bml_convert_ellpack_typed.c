#include "../../macros.h"
#include "../../typed.h"
#include "../bml_getters.h"
#include "../bml_introspection.h"
#include "../bml_logger.h"
#include "../bml_setters.h"
#include "bml_allocate_ellpack.h"
#include "bml_types_ellpack.h"

#include <complex.h>

bml_matrix_ellpack_t *TYPED_FUNC(
    bml_convert_ellpack) (
    bml_matrix_t * A,
    bml_matrix_precision_t matrix_precision,
    int M,
    bml_distribution_mode_t distrib_mode)
{
    int N = bml_get_N(A);

    bml_matrix_type_t A_type = bml_get_type(A);
    if (A_type == ellpack) {
        bml_matrix_ellpack_t *A_ellpack = (bml_matrix_ellpack_t*)A;
        int *A_nnz = A_ellpack->nnz;
        int *A_index = A_ellpack->index;
        REAL_T *A_value = A_ellpack->value;

#pragma omp target update from(A_nnz[:N], A_index[:N*M], A_value[:N*M])
    }

    if (N < 0)
    {
        LOG_ERROR("A is not intialized\n");
    }

    bml_matrix_ellpack_t *B =
        bml_zero_matrix_ellpack(matrix_precision, N, M, distrib_mode);

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            bml_set_element(B, i, j, bml_get(A, i, j));
        }
    }

    int *B_nnz = B->nnz;
    int *B_index = B->index;
    REAL_T *B_value = B->value;
#pragma omp target update to(B_nnz[:N], B_index[:N*M], B_value[:N*M])
    return B;
}
