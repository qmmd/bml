#include "../../macros.h"
#include "../../typed.h"
#include "../bml_parallel.h"
#include "../bml_types.h"
#include "bml_parallel_dense.h"
#include "bml_types_dense.h"
#include "bml_allocate_dense.h"
#include "../bml_logger.h"

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
#if defined(DO_MPI) && defined(BML_MPI_NONDIST) 
    int myRank = bml_getMyRank();

    int N = A->N;

    REAL_T *A_matrix = A->matrix;

    bml_domain_t *A_domain = A->domain;

    MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                   A_matrix, A_domain->localElements,
                   A_domain->localDispl, REAL_MPI_TYPE, ccomm);
#endif
}

#ifdef DO_MPI

void TYPED_FUNC(
    bml_mpi_type_create_struct_dense) (
    bml_matrix_dense_t * A,
    MPI_Datatype * newtype)
{
    MPI_Aint baseaddr;
    MPI_Aint addr0, addr1, addr2;
    MPI_Get_address(A, &baseaddr);
    MPI_Get_address(A->matrix, &addr0);

    MPI_Datatype dtype[1];
    dtype[0] = MPI_T;

    int blength[1];
    blength[0] = A->N * A->N;

    MPI_Aint displ[0];
    displ[0] = addr0 - baseaddr;
    int mpiret = MPI_Type_create_struct(1, blength, displ, dtype, newtype);
    if (mpiret != MPI_SUCCESS)
        LOG_ERROR("MPI_Type_create_struct failed!");
    mpiret = MPI_Type_commit(newtype);
    if (mpiret != MPI_SUCCESS)
        LOG_ERROR("MPI_Type_commit failed!");
}

void TYPED_FUNC(
    bml_mpi_send_dense) (
    bml_matrix_dense_t * A,
    const int dst,
    MPI_Comm comm)
{
    MPI_Send(A->matrix, A->N * A->N, MPI_T, dst, 222, comm);
}

void TYPED_FUNC(
    bml_mpi_recv_dense) (
    bml_matrix_dense_t * A,
    const int src,
    MPI_Comm comm)
{
    MPI_Status status;
    MPI_Recv(A->matrix, A->N * A->N, MPI_T, src, 222, comm, &status);
}

bml_matrix_dense_t
    * TYPED_FUNC(bml_mpi_recv_matrix_dense) (int N, int M,
                                             const int src, MPI_Comm comm)
{
    bml_matrix_dimension_t matrix_dimension = { N, N, M };
    bml_matrix_dense_t *A_bml =
        TYPED_FUNC(bml_zero_matrix_dense) (matrix_dimension, sequential);

    bml_mpi_recv_dense(A_bml, src, comm);

    return A_bml;
}

#endif
