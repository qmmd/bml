#include "../../macros.h"
#include "../../typed.h"
#include "../lapack.h"
#include "../bml_allocate.h"
#include "../bml_logger.h"
#include "../bml_types.h"
#include "../bml_export.h"
#include "bml_allocate_distributed2d.h"
#include "bml_export_distributed2d.h"
#include "bml_types_distributed2d.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

/** Convert a distributed bml matrix into a dense matrix
 *  on MPI task 0
 *
 * \ingroup convert_group
 *
 * \param A The bml matrix
 * \return The dense matrix
 */
void *TYPED_FUNC(
    bml_export_to_dense_distributed2d) (
    bml_matrix_distributed2d_t * A_bml,
    bml_dense_order_t order)
{
    int myrank = A_bml->mpitask;
    int ntasks = A_bml->ntasks;

    int N = A_bml->N;
    int n = A_bml->n;

    REAL_T *A_dense = NULL;
    if (myrank == 0)
        A_dense = bml_allocate_memory(sizeof(REAL_T) * N * N);

    // export local submatrix
    REAL_T *sendbuf = bml_export_to_dense(A_bml->matrix, order);
    // gather data on MPI task 0
    int tag = 0;
    MPI_Status status;
    if (myrank != 0)
    {
        int dest = 0;
        MPI_Send(sendbuf, n * n, MPI_T, dest, tag, A_bml->comm);
    }
    else
    {
        // copy local data into A_dense
        C_BLAS(LACPY) ("A", &n, &n, sendbuf, &n, A_dense, &N);
        REAL_T **recvbuf =
            bml_allocate_memory((ntasks - 1) * sizeof(REAL_T *));
        MPI_Request *request =
            bml_allocate_memory((ntasks - 1) * sizeof(MPI_Request));
        for (int src = 1; src < ntasks; src++)
        {
            recvbuf[src - 1] = bml_allocate_memory(n * n * sizeof(REAL_T));
            MPI_Irecv(recvbuf[src - 1], n * n, MPI_T, src, tag, A_bml->comm,
                      &request[src - 1]);
        }
        int coords[2];
        for (int src = 1; src < ntasks; src++)
        {
            MPI_Cart_coords(A_bml->comm, src, 2, coords);
            int ip = coords[0];
            int jp = coords[1];
            MPI_Wait(&request[src - 1], &status);
            C_BLAS(LACPY) ("A", &n, &n, recvbuf[src - 1], &n,
                           A_dense + n * N * ip + n * jp, &N);
        }
        for (int src = 1; src < ntasks; src++)
        {
            bml_free_memory(recvbuf[src - 1]);
        }
        bml_free_memory(request);
        bml_free_memory(recvbuf);
    }

    return A_dense;
}
