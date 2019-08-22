#include "../bml_logger.h"
#include "../bml_parallel.h"
#include "../bml_types.h"
#include "bml_parallel_ellsort.h"
#include "bml_types_ellsort.h"

#include <stdlib.h>
#include <string.h>

/** Gather bml matrix across MPI ranks.
 *
 *  \ingroup parallel_group
 *
 *  \param A The matrix
 */
void
bml_allGatherVParallel_ellsort(
    bml_matrix_ellsort_t * A)
{

    switch (A->matrix_precision)
    {
        case single_real:
            bml_allGatherVParallel_ellsort_single_real(A);
            break;
        case double_real:
            bml_allGatherVParallel_ellsort_double_real(A);
            break;
        case single_complex:
            bml_allGatherVParallel_ellsort_single_complex(A);
            break;
        case double_complex:
            bml_allGatherVParallel_ellsort_double_complex(A);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
}
