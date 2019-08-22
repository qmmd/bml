#include "../bml_logger.h"
#include "../bml_parallel.h"
#include "../bml_types.h"
#include "bml_parallel_dense.h"
#include "bml_types_dense.h"

#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/** Gather pieces of matrix across MPI ranks.
 *
 *  \ingroup parallel_group
 *
 *  \param A The matrix
 */
void
bml_allGatherVParallel_dense(
    bml_matrix_dense_t * A)
{
    switch (A->matrix_precision)
    {
        case single_real:
            bml_allGatherVParallel_dense_single_real(A);
            break;
        case double_real:
            bml_allGatherVParallel_dense_double_real(A);
            break;
        case single_complex:
            bml_allGatherVParallel_dense_single_complex(A);
            break;
        case double_complex:
            bml_allGatherVParallel_dense_double_complex(A);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
}
