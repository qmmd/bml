#include "../../macros.h"
#include "../../typed.h"
#include "../bml_parallel.h"
#include "../bml_types.h"
#include "bml_parallel_ellsort.h"
#include "bml_types_ellsort.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef DO_MPI
#include <mpi.h>
#endif

/** Gather a bml matrix across MPI ranks.
 *
 *  \ingroup parallel_group
 *
 *  \param A The matrix A
 */
void TYPED_FUNC(
    bml_allGatherVParallel_ellsort) (
    bml_matrix_ellsort_t * A)
{
    int myRank = bml_getMyRank();
    int nRanks = bml_getNRanks();

    int N = A->N;
    int M = A->M;

    int *A_nnz = (int *) A->nnz;
    int *A_index = (int *) A->index;
    bml_domain_t *A_domain = (bml_domain_t *) A->domain;

    REAL_T *A_value = (REAL_T *) A->value;

#ifdef DO_MPI
/*
    for (int i = 0; i < nRanks; i++)
    {
      printf("allgatherv %d: rank %d localRow %d %d %d\n", myRank, i, A_domain->localRowMin[i],
        A_domain->localRowMax[i], A_domain->localRowExtent[i]);
      printf("allgatherv %d before: rank %d nnz  %d %d %d %d\n", myRank, i, A_nnz[0],
        A_nnz[3071], A_nnz[3072], A_nnz[6143]);
    }
*/

    // Number of non-zeros per row
    MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                   A_nnz, A_domain->localRowExtent,
                   A_domain->localRowMin, MPI_INT, ccomm);

    // Indeces
    MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                   A_index, A_domain->localElements, A_domain->localDispl,
                   MPI_INT, ccomm);

    // Values
    MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                   A_value, A_domain->localElements, A_domain->localDispl,
                   REAL_MPI_TYPE, ccomm);

/*
    for (int i = 0; i < nRanks; i++)
    {
      printf("allgatherv %d: rank %d localRow %d %d %d\n", myRank, i, A_domain->localRowMin[i],
        A_domain->localRowMax[i], A_domain->localRowExtent[i]);
      printf("allgatherv %d after: rank %d nnz  %d %d %d %d\n", myRank, i, A_nnz[0],
        A_nnz[3071], A_nnz[3072], A_nnz[6143]);
    }
*/
#endif

}
