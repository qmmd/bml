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

#ifdef DO_MPI
void
bml_mpi_type_create_struct_dense(
    bml_matrix_dense_t * A,
    MPI_Datatype * newtype)
{

    switch (A->matrix_precision)
    {
        case single_real:
            bml_mpi_type_create_struct_dense_single_real(A, newtype);
            break;
        case double_real:
            bml_mpi_type_create_struct_dense_double_real(A, newtype);
            break;
        case single_complex:
            bml_mpi_type_create_struct_dense_single_complex(A, newtype);
            break;
        case double_complex:
            bml_mpi_type_create_struct_dense_double_complex(A, newtype);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
}

void
bml_mpi_send_dense(
    bml_matrix_dense_t * A,
    const int dst,
    MPI_Comm comm)
{
    switch (A->matrix_precision)
    {
        case single_real:
            bml_mpi_send_dense_single_real(A, dst, comm);
            break;
        case double_real:
            bml_mpi_send_dense_double_real(A, dst, comm);
            break;
        case single_complex:
            bml_mpi_send_dense_single_complex(A, dst, comm);
            break;
        case double_complex:
            bml_mpi_send_dense_double_complex(A, dst, comm);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
}

void
bml_mpi_recv_dense(
    bml_matrix_dense_t * A,
    const int dst,
    MPI_Comm comm)
{
    switch (A->matrix_precision)
    {
        case single_real:
            bml_mpi_recv_dense_single_real(A, dst, comm);
            break;
        case double_real:
            bml_mpi_recv_dense_double_real(A, dst, comm);
            break;
        case single_complex:
            bml_mpi_recv_dense_single_complex(A, dst, comm);
            break;
        case double_complex:
            bml_mpi_recv_dense_double_complex(A, dst, comm);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
}

bml_matrix_dense_t *
bml_mpi_recv_matrix_dense(
    bml_matrix_precision_t matrix_precision,
    int N,
    int M,
    const int src,
    MPI_Comm comm)
{
    switch (matrix_precision)
    {
        case single_real:
            return bml_mpi_recv_matrix_dense_single_real(N, M, src, comm);
            break;
        case double_real:
            return bml_mpi_recv_matrix_dense_double_real(N, M, src, comm);
            break;
        case single_complex:
            return bml_mpi_recv_matrix_dense_single_complex(N, M, src, comm);
            break;
        case double_complex:
            return bml_mpi_recv_matrix_dense_double_complex(N, M, src, comm);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
}

void
bml_mpi_bcast_matrix_dense(
    bml_matrix_dense_t * A,
    const int root,
    MPI_Comm comm)
{
    // create MPI data type to avoid multiple messages
    MPI_Datatype mpi_data_type;
    bml_mpi_type_create_struct_dense(A, &mpi_data_type);

    MPI_Bcast(A, 1, mpi_data_type, root, comm);

    MPI_Type_free(&mpi_data_type);
}

#endif
