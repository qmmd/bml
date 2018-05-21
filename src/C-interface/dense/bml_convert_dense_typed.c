#include "../../macros.h"
#include "../../typed.h"
#include "bml_allocate_dense.h"
#include "bml_getters.h"
#include "bml_introspection.h"
#include "bml_logger.h"
#include "bml_types_dense.h"

#include <complex.h>

bml_matrix_dense_t *TYPED_FUNC(
    bml_convert_dense) (
    const bml_matrix_t * A,
    const bml_matrix_precision_t matrix_precision,
    const bml_distribution_mode_t distrib_mode)
{
    int N = bml_get_N(A);

    if (N < 0)
    {
        LOG_ERROR("A is not initialized\n");
    }

    bml_matrix_dense_t *B =
        bml_zero_matrix_dense(matrix_precision, N, distrib_mode);
    REAL_T *Bij = (REAL_T *) B->matrix;

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            Bij[LINEARINDEX(i, j, N, N)] = *(REAL_T *) bml_get(A, i, j);
        }
    }

    return B;
}
