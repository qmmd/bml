#ifdef BML_USE_MAGMA
#include "magma_v2.h"
#else
#include "../blas.h"
#include "../lapack.h"
#endif
#include "../../typed.h"
#include "../bml_allocate.h"
#include "../bml_copy.h"
#include "../bml_inverse.h"
#include "../bml_logger.h"
#include "../bml_parallel.h"
#include "../bml_types.h"
#include "bml_copy_dense.h"
#include "bml_inverse_dense.h"
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
    bml_matrix_dense_t * A)
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
    int lda = A->ld;
    int *ipiv = bml_allocate_memory(N * sizeof(int));
#ifdef BML_USE_MAGMA
    MAGMAGPU(getrf) (M, N, A->matrix, A->ld, ipiv, &info);
    if (info != 0)
        LOG_ERROR("ERROR in getrf_gpu");
    int nb = magma_get_dgetri_nb(N);
    int lwork = N * nb;
    MAGMA_T *dwork;
    MAGMA(malloc) (&dwork, lwork * sizeof(MAGMA_T));
    MAGMAGPU(getri) (N, A->matrix, A->ld, ipiv, dwork, lwork, &info);
    if (info != 0)
        LOG_ERROR("ERROR in getri_gpu");
    magma_free(dwork);
#else
    int lwork = N * N;
    REAL_T *work = bml_allocate_memory(lwork * sizeof(REAL_T));

#ifdef NOBLAS
    LOG_ERROR("No BLAS library");
#else
    C_BLAS(GETRF) (&M, &N, A->matrix, &lda, ipiv, &info);
    C_BLAS(GETRI) (&N, A->matrix, &N, ipiv, work, &lwork, &info);
#endif
    bml_free_memory(work);
#endif /* BML_USE_MAGMA */
    bml_free_memory(ipiv);
}
