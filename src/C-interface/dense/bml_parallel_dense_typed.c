#include "../../macros.h"
#include "../../typed.h"
#include "bml_parallel.h"
#include "bml_parallel_dense.h"
#include "bml_types.h"
#include "bml_types_dense.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/** Gather a bml matrix across MPI ranks.
 *
 *  \ingroup parallel_group
 *
 *  \param A The matrix
 */
void TYPED_FUNC(
    bml_allGatherVParallel_dense) (
    bml_matrix_dense_t * A)
{
#ifdef DO_MPI
    int myRank = bml_getMyRank();

    int N = A->N;

    REAL_T *A_matrix = A->matrix;

    bml_domain_t *A_domain = A->domain;

    MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                   A_matrix, A_domain->localElements,
                   A_domain->localDispl, REAL_MPI_TYPE, ccomm);
#endif
}
