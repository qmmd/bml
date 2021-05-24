#include "../../macros.h"
#include "../../typed.h"
#include "../bml_parallel.h"
#include "../bml_types.h"
#include "bml_parallel_ellpack.h"
#include "bml_types_ellpack.h"
#include "bml_allocate_ellpack.h"
#include "../bml_logger.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#ifdef _OPENMP
#include <omp.h>
#endif

//#ifdef DO_MPI
#include <mpi.h>
//#endif

/** Gather a bml matrix across MPI ranks.
 *
 *  \ingroup parallel_group
 *
 *  \param A The matrix A
 */
void TYPED_FUNC(
    bml_allGatherVParallel_ellpack) (
    bml_matrix_ellpack_t * A)
{
    int myRank = bml_getMyRank();
    int nRanks = bml_getNRanks();

    int N = A->N;
    int M = A->M;

    int *A_nnz = (int *) A->nnz;
    int *A_index = (int *) A->index;
    bml_domain_t *A_domain = (bml_domain_t *) A->domain;

    REAL_T *A_value = (REAL_T *) A->value;

#if defined(DO_MPI) && defined(BML_MPI_NONDIST) 
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

#ifdef DO_MPI

void TYPED_FUNC(
    bml_mpi_type_create_struct_ellpack) (
    bml_matrix_ellpack_t * A,
    MPI_Datatype * newtype)
{
    MPI_Aint baseaddr;
    MPI_Aint addr0, addr1, addr2;
    MPI_Get_address(A, &baseaddr);
    MPI_Get_address(A->value, &addr0);
    MPI_Get_address(A->index, &addr1);
    MPI_Get_address(A->nnz, &addr2);

    MPI_Datatype dtype[3];
    dtype[0] = MPI_T;
    dtype[1] = MPI_INT;
    dtype[2] = MPI_INT;

    int blength[3];
    blength[0] = A->N * A->M;
    blength[1] = A->N * A->M;
    blength[2] = A->N;

    MPI_Aint displ[3];
    displ[0] = addr0 - baseaddr;
    displ[1] = addr1 - baseaddr;
    displ[2] = addr2 - baseaddr;
    int mpiret = MPI_Type_create_struct(3, blength, displ, dtype, newtype);
    if (mpiret != MPI_SUCCESS)
        LOG_ERROR("MPI_Type_create_struct failed!");
    mpiret = MPI_Type_commit(newtype);
    if (mpiret != MPI_SUCCESS)
        LOG_ERROR("MPI_Type_commit failed!");
}

void TYPED_FUNC(
    bml_mpi_send_ellpack) (
    bml_matrix_ellpack_t * A,
    const int dst,
    MPI_Comm comm)
{
    // create MPI data type to avoid multiple messages
    MPI_Datatype mpi_data_type;
    bml_mpi_type_create_struct_ellpack(A, &mpi_data_type);

    MPI_Send(A, 1, mpi_data_type, dst, 111, comm);

    MPI_Type_free(&mpi_data_type);
}

void TYPED_FUNC(
    bml_mpi_recv_ellpack) (
    bml_matrix_ellpack_t * A,
    const int src,
    MPI_Comm comm)
{
    // create MPI data type to avoid multiple messages
    MPI_Datatype mpi_data_type;
    bml_mpi_type_create_struct_ellpack(A, &mpi_data_type);

    MPI_Status status;
    MPI_Recv(A, 1, mpi_data_type, src, 111, comm, &status);

    MPI_Type_free(&mpi_data_type);
}

/*
 * Return BML matrix from data received from MPI task src
 */
bml_matrix_ellpack_t
    * TYPED_FUNC(bml_mpi_recv_matrix_ellpack) (int N, int M,
                                               const int src, MPI_Comm comm)
{
    bml_matrix_ellpack_t *A_bml =
        TYPED_FUNC(bml_zero_matrix_ellpack) (N, M, sequential);

    bml_mpi_recv_ellpack(A_bml, src, comm);

    return A_bml;
}

#endif
