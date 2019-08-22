#include "../../macros.h"
#include "../../typed.h"
#include "../bml_getters.h"
#include "../bml_introspection.h"
#include "../bml_logger.h"
#include "../bml_setters.h"
#include "bml_allocate_ellsort.h"
#include "bml_types_ellsort.h"

#include <complex.h>

bml_matrix_ellsort_t *TYPED_FUNC(
    bml_convert_ellsort) (
    bml_matrix_t * A,
    bml_matrix_precision_t matrix_precision,
    int M,
    bml_distribution_mode_t distrib_mode)
{
    int N = bml_get_N(A);

    if (N < 0)
    {
        LOG_ERROR("A is not intialized\n");
    }

    bml_matrix_ellsort_t *B =
        bml_zero_matrix_ellsort(matrix_precision, N, M, distrib_mode);

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            bml_set_element(B, i, j, bml_get(A, i, j));
        }
    }

    return B;
}
