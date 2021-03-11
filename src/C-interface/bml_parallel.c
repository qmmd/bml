#include "bml_parallel.h"
#include "bml_introspection.h"
#include "bml_logger.h"
#include "dense/bml_parallel_dense.h"
#include "ellpack/bml_parallel_ellpack.h"
#include "ellsort/bml_parallel_ellsort.h"
#include "ellblock/bml_parallel_ellblock.h"
#include "csr/bml_parallel_csr.h"
#ifdef DO_MPI
#include "distributed2d/bml_allocate_distributed2d.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <assert.h>

#ifdef DO_MPI
#include <mpi.h>
#endif

static int myRank = 0;
static int nRanks = 1;
#ifdef DO_MPI
static MPI_Request *requestList;
MPI_Comm ccomm;
static int *rUsed;
#endif

/** Initialize.
 *
 * \ingroup init_group_C
 *
 * \param argc Number of args
 * \param argv Args
 */

/** Get number of MPI ranks.
 */
int
bml_getNRanks(
    void)
{
    return nRanks;
}

/** Get local MPI rank.
 */
int
bml_getMyRank(
    void)
{
    return myRank;
}

/* \details
 *
 * For now this is just a check for rank 0 but in principle it could be more
 * complex.  It is also possible to suppress practically all output by causing
 * this function to return 0 for all ranks.
 */
int
bml_printRank(
    void)
{
    if (myRank == 0)
        return 1;
    return 0;
}

#ifdef DO_MPI
void
bml_initParallel(
    MPI_Comm comm)
{
    ccomm = comm;
    MPI_Comm_rank(ccomm, &myRank);
    MPI_Comm_size(ccomm, &nRanks);

    if (bml_printRank())
        printf("MPI started in bml with %d ranks\n", nRanks);

    bml_setcomm_distributed2d(comm);

    requestList = (MPI_Request *) malloc(nRanks * sizeof(MPI_Request));
    rUsed = (int *) malloc(nRanks * sizeof(int));
    for (int i = 0; i < nRanks; i++)
    {
        rUsed[i] = 0;
    }
}
#endif

void
bml_initParallelF(
    int fcomm)
{
#ifdef DO_MPI
    MPI_Comm comm = MPI_Comm_f2c(fcomm);
    bml_initParallel(comm);
#endif
}

void
bml_shutdownParallelF(
    void)
{
    bml_shutdownParallel();
}

void
bml_shutdownParallel(
    void)
{
#ifdef DO_MPI
    free(requestList);
    free(rUsed);
#endif
}

void
bml_barrierParallel(
    void)
{
#ifdef DO_MPI
    MPI_Barrier(ccomm);
#endif
}

void
bml_sumRealReduce(
    double *value)
{
#ifdef DO_MPI
    double sLocal[1], sGlobal[1];

    sLocal[0] = *value;
    MPI_Allreduce(sLocal, sGlobal, 1, MPI_DOUBLE, MPI_SUM, ccomm);
    *value = sGlobal[0];
#endif
}

void
bml_minRealReduce(
    double *value)
{
#ifdef DO_MPI
    double sLocal[1], sGlobal[1];

    sLocal[0] = *value;
    MPI_Allreduce(sLocal, sGlobal, 1, MPI_DOUBLE, MPI_MIN, ccomm);
    *value = sGlobal[0];
#endif
}

void
bml_maxRealReduce(
    double *value)
{
#ifdef DO_MPI
    double sLocal[1], sGlobal[1];

    sLocal[0] = *value;
    MPI_Allreduce(sLocal, sGlobal, 1, MPI_DOUBLE, MPI_MAX, ccomm);
    *value = sGlobal[0];
#endif
}

/** Exchange pieces of matrix across MPI ranks.
 *
 * \param A Matrix A
 */
void
bml_allGatherVParallel(
    bml_matrix_t * A)
{
    switch (bml_get_type(A))
    {
        case dense:
            bml_allGatherVParallel_dense(A);
            break;
        case ellpack:
            bml_allGatherVParallel_ellpack(A);
            break;
        case ellsort:
            bml_allGatherVParallel_ellsort(A);
            break;
        case ellblock:
            LOG_ERROR
                ("bml_allGatherVParallel function not implemented for ellblock\n");
        default:
            LOG_ERROR("unknown matrix type\n");
            break;
    }
}

#ifdef DO_MPI
void
bml_mpi_send(
    bml_matrix_t * A,
    const int dst,
    MPI_Comm comm)
{
    switch (bml_get_type(A))
    {
        case dense:
            bml_mpi_send_dense(A, dst, comm);
            break;
        case ellpack:
            bml_mpi_send_ellpack(A, dst, comm);
            break;
        case ellsort:
            bml_mpi_send_ellsort(A, dst, comm);
            break;
        case ellblock:
            bml_mpi_send_ellblock(A, dst, comm);
            break;
        case csr:
            bml_mpi_send_csr(A, dst, comm);
            break;
        default:
            LOG_ERROR("unknown matrix type\n");
            break;
    }
}

void
bml_mpi_recv(
    bml_matrix_t * A,
    const int src,
    MPI_Comm comm)
{
    switch (bml_get_type(A))
    {
        case dense:
            bml_mpi_recv_dense(A, src, comm);
            break;
        case ellpack:
            bml_mpi_recv_ellpack(A, src, comm);
            break;
        case ellsort:
            bml_mpi_recv_ellsort(A, src, comm);
            break;
        case ellblock:
            bml_mpi_recv_ellblock(A, src, comm);
            break;
        case csr:
            bml_mpi_recv_csr(A, src, comm);
            break;
        default:
            LOG_ERROR("unknown matrix type\n");
            break;
    }
}

bml_matrix_t *
bml_mpi_recv_matrix(
    bml_matrix_type_t matrix_type,
    bml_matrix_precision_t matrix_precision,
    int N,
    int M,
    const int src,
    MPI_Comm comm)
{
    switch (matrix_type)
    {
        case dense:
            return bml_mpi_recv_matrix_dense(matrix_precision, N, M, src,
                                             comm);
            break;
        case ellpack:
            return bml_mpi_recv_matrix_ellpack(matrix_precision, N, M, src,
                                               comm);
            break;
        case ellsort:
            return bml_mpi_recv_matrix_ellsort(matrix_precision, N, M, src,
                                               comm);
            break;
        case ellblock:
            return bml_mpi_recv_matrix_ellblock(matrix_precision, N, M, src,
                                                comm);
            break;
        case csr:
            return bml_mpi_recv_matrix_csr(matrix_precision, N, M, src, comm);
            break;
        default:
            LOG_ERROR("unknown matrix type\n");
            break;
    }
}

#endif
