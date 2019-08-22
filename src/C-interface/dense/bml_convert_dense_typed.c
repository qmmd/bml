#ifdef BML_USE_MAGMA
#include "magma_v2.h"
#endif

#include "../../macros.h"
#include "../../typed.h"
#include "../bml_getters.h"
#include "../bml_introspection.h"
#include "../bml_logger.h"
#include "bml_allocate_dense.h"
#include "bml_types_dense.h"

#include <complex.h>

bml_matrix_dense_t *TYPED_FUNC(
    bml_convert_dense) (
    bml_matrix_t * A,
    bml_matrix_precision_t matrix_precision,
    bml_distribution_mode_t distrib_mode)
{
    int N = bml_get_N(A);

    if (N < 0)
    {
        LOG_ERROR("A is not initialized\n");
    }

    bml_matrix_dimension_t matrix_dimension = { N, N, N };
    bml_matrix_dense_t *B =
        bml_zero_matrix_dense(matrix_precision, matrix_dimension,
                              distrib_mode);
#ifdef BML_USE_MAGMA
    REAL_T *Bij = calloc(N * N, sizeof(REAL_T));
#else
    REAL_T *Bij = (REAL_T *) B->matrix;
#endif

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            Bij[ROWMAJOR(i, j, N, N)] = *(REAL_T *) bml_get(A, i, j);
        }
    }
#ifdef BML_USE_MAGMA
    MAGMA(setmatrix) (N, N, (MAGMA_T *) Bij, N, B->matrix, B->ld, B->queue);
    free(Bij);
#endif
    return B;
}
