#include "../../macros.h"
#include "../../typed.h"
#include "../bml_multiply.h"
#include "../bml_allocate.h"
#include "../bml_logger.h"
#include "../bml_copy.h"
#include "../bml_parallel.h"
#include "bml_allocate_distributed2d.h"
#include "bml_import_distributed2d.h"
#include "bml_types_distributed2d.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <mpi.h>

/** Matrix multiply.
 *
 * C = alpha * A * B + beta * C
 *
 *  \ingroup multiply_group
 *
 *  \param A Matrix A
 *  \param B Matrix B
 *  \param C Matrix C
 *  \param alpha Scalar factor multiplied by A * B
 *  \param beta Scalar factor multiplied by C
 *  \param threshold Used for sparse multiply
 */
void TYPED_FUNC(
    bml_multiply_distributed2d) (
    bml_matrix_distributed2d_t * A,
    bml_matrix_distributed2d_t * B,
    bml_matrix_distributed2d_t * C,
    double alpha,
    double beta,
    double threshold)
{
    // make a copy of A and B local submatrices
    bml_matrix_t *Atmp1 = NULL;
    bml_matrix_t *Atmp2 = NULL;
    Atmp1 = bml_copy_new(A->matrix);
    Atmp2 = bml_copy_new(A->matrix);

    bml_matrix_t *Btmp1 = NULL;
    bml_matrix_t *Btmp2 = NULL;
    Btmp1 = bml_copy_new(B->matrix);
    Btmp2 = bml_copy_new(B->matrix);

    // shift all submatrices A(i,j) to the left by i steps
    int src, dst;
    MPI_Cart_shift(A->comm, 1, A->myprow, &src, &dst);
    bml_mpi_send(Atmp1, dst, A->comm);
    bml_mpi_recv(Atmp2, src, A->comm);

    // shift all submatrices B(i,j) up by j steps
    MPI_Cart_shift(B->comm, 0, B->mypcol, &src, &dst);
    bml_mpi_send(Btmp1, dst, B->comm);
    bml_mpi_recv(Btmp2, src, B->comm);

    // perform local submatrices multiplication
    bml_multiply(Atmp2, Btmp2, C->matrix, alpha, beta, threshold);

    for (int k = 1; k < C->npcols; k++)
    {
        //swap pointers
        bml_matrix_t *tmp = Atmp1;
        Atmp1 = Atmp2;
        Atmp2 = tmp;
        tmp = Btmp1;
        Btmp1 = Btmp2;
        Btmp2 = tmp;
        // Move submatrices A(i,j) one step to left
        MPI_Cart_shift(C->comm, 1, 1, &src, &dst);
        bml_mpi_send(Atmp1, dst, A->comm);
        bml_mpi_recv(Atmp2, src, A->comm);

        // Move submatrices B(i,j) one step up
        MPI_Cart_shift(C->comm, 0, 1, &src, &dst);
        bml_mpi_send(Btmp1, dst, B->comm);
        bml_mpi_recv(Btmp2, src, B->comm);

        // perform local submatrices multiplication
        bml_multiply(Atmp2, Btmp2, C->matrix, alpha, 1., threshold);
    }

//    bml_deallocate(&Atmp1);
//    bml_deallocate(&Atmp2);
//    bml_deallocate(&Btmp1);
//    bml_deallocate(&Btmp2);
}
