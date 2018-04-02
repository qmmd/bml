#include "../blas.h"
#include "../lapack.h"
#include "../../typed.h"
#include "../bml_logger.h"
#include "bml_copy.h"
#include "bml_inverse.h"
#include "bml_allocate.h"
#include "bml_parallel.h"
#include "bml_copy_dense.h"
#include "bml_inverse_dense.h"
#include "bml_types.h"
#include "bml_types_dense.h"

#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/** Matrix inverse.
 *
 * \f$ B \leftarrow A^-1 f$
 *
 * \ingroup inverse_group
 *
 * \param A Matrix A
 * \param B Inverse of Matrix A
 */
bml_matrix_dense_t *TYPED_FUNC(
    bml_inverse_dense) (
    const bml_matrix_dense_t * A)
{
    bml_matrix_dense_t *B = TYPED_FUNC(bml_copy_dense_new) (A);

    TYPED_FUNC(bml_inverse_inplace_dense) (B);

    return B;
}

/** Matrix inverse inplace.
 *
 * \f$ A \leftarrow A^-1 f$
 *
 * \ingroup inverse_group
 *
 * \param A Matrix A
 */
void TYPED_FUNC(
    bml_inverse_inplace_dense) (
    bml_matrix_dense_t * A)
{
    int info;
    int M = A->N;
    int N = A->N;
    int lda = A->N;
    int lwork = N * N;
    int *ipiv = bml_allocate_memory(N * sizeof(int));
    REAL_T *work = bml_allocate_memory(lwork * sizeof(REAL_T));

#ifdef NOBLAS
    LOG_ERROR("No BLAS library");
#else
    C_BLAS(GETRF) (&M, &N, A->matrix, &lda, ipiv, &info);
    C_BLAS(GETRI) (&N, A->matrix, &N, ipiv, work, &lwork, &info);
#endif

    bml_free_memory(ipiv);
    bml_free_memory(work);
}
