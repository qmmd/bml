#include "../../macros.h"
#include "../../typed.h"
#include "../lapack.h"
#include "../bml_allocate.h"
#include "../bml_logger.h"
#include "../bml_import.h"
#include "bml_allocate_distributed2d.h"
#include "bml_import_distributed2d.h"
#include "bml_types_distributed2d.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

/** Convert and distribute a dense matrix on MPI task 0
 *  into a bml block matrix on each task.
 *
 * \ingroup convert_group
 *
 * \param N The number of rows/columns
 * \param matrix_precision The real precision
 * \param A The dense matrix
 * \param threshold The matrix element magnited threshold
 * \return The bml matrix
 */
bml_matrix_distributed2d_t *TYPED_FUNC(
    bml_import_from_dense_distributed2d) (
    bml_matrix_type_t matrix_type,
    bml_dense_order_t order,
    int N,
    void *A,
    double threshold,
    int M)
{
    bml_matrix_distributed2d_t *A_bml =
        bml_allocate_memory(sizeof(bml_matrix_distributed2d_t));
    // setup A_bml struct
    A_bml->matrix_precision = MATRIX_PRECISION;
    bml_setup_distributed2d(N, A_bml);
    assert(A_bml->comm != MPI_COMM_NULL);

    // local submatrix dimensions
    int n = A_bml->n;
    int myrank = A_bml->mpitask;
    int ntasks = A_bml->ntasks;
    int m = M / (int) sqrt(ntasks);
    assert(m <= n);

    REAL_T *recvbuf = bml_allocate_memory(n * n * sizeof(REAL_T));
    int tag = 0;
    MPI_Status status;
    if (myrank == 0)            // send data to other MPI tasks
    {
        int coords[2];
        REAL_T **sendbuf =
            bml_allocate_memory((ntasks - 1) * sizeof(REAL_T *));
        MPI_Request *request =
            bml_allocate_memory((ntasks - 1) * sizeof(MPI_Request));
        for (int dest = 1; dest < ntasks; dest++)
        {
            // find out which block of data to sent
            MPI_Cart_coords(A_bml->comm, dest, 2, coords);
            int ip = coords[0];
            int jp = coords[1];
            // pack data into buffer
            sendbuf[dest - 1] = bml_allocate_memory(n * n * sizeof(REAL_T));
            REAL_T *array = (REAL_T *) (A);
            C_BLAS(LACPY) ("A", &n, &n, array + N * n * ip + n * jp, &N,
                           sendbuf[dest - 1], &n);
            MPI_Isend(sendbuf[dest - 1], n * n, MPI_T, dest, tag, A_bml->comm,
                      &request[dest - 1]);
        }
        // put local data into recvbuf
        C_BLAS(LACPY) ("A", &n, &n, A, &N, recvbuf, &n);
        for (int dest = 1; dest < ntasks; dest++)
        {
            MPI_Wait(&request[dest - 1], &status);
            bml_free_memory(sendbuf[dest - 1]);
        }

        bml_free_memory(request);
        bml_free_memory(sendbuf);
    }
    else                        // receive data from task 0
    {
        int source = 0;
        MPI_Recv(recvbuf, n * n, MPI_T, source, tag, A_bml->comm,
                 MPI_STATUS_IGNORE);
    }

    // build local submatrix with received data
    A_bml->matrix =
        bml_import_from_dense(matrix_type, MATRIX_PRECISION, order, n, m,
                              recvbuf, threshold, sequential);

    bml_free_memory(recvbuf);

    return A_bml;
}
