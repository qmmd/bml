#include "../../macros.h"
#include "../../typed.h"
#include "../bml_getters.h"
#include "../bml_introspection.h"
#include "../bml_logger.h"
#include "../bml_setters.h"
#include "../bml_allocate.h"
#include "bml_allocate_csr.h"
#include "bml_types_csr.h"

#include <complex.h>

/** Convert from dense matrix to a csr matrix.
 *
 *  \ingroup convert
 *
 *  \param A The dense matrix A
 *  \return The csr format of A
 */
bml_matrix_csr_t *TYPED_FUNC(
    bml_convert_csr) (
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

    bml_matrix_csr_t *B =
        bml_zero_matrix_csr(matrix_precision, N, M, distrib_mode);

    for (int i = 0; i < N; i++)
    {
        REAL_T *row = bml_get_row(A, i);
        bml_set_row(B, i, row, 0.0);
        bml_free_memory(row);
    }

    return B;
}
