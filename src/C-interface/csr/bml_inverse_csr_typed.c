#include "bml_allocate_csr.h"
#include "../bml_allocate.h"
#include "bml_copy_csr.h"
#include "bml_export_csr.h"
#include "bml_inverse_csr.h"
#include "../bml_inverse.h"
#include "bml_types_csr.h"
#include "../bml_types.h"
#include "../bml_utilities.h"
#include "../dense/bml_allocate_dense.h"
#include "../dense/bml_import_dense.h"
#include "../dense/bml_inverse_dense.h"
#include "../dense/bml_types_dense.h"
#include "../csr/bml_import_csr.h"
#include "../../macros.h"
#include "../../typed.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/** Matrix inverse.
 *
 *  \ingroup inverse_group
 *
 *  \param A The matrix A
 *  \return The inverse of matrix A
 */
bml_matrix_csr_t *TYPED_FUNC(
    bml_inverse_csr) (
    bml_matrix_csr_t * A)
{
    const int N = A->N_;
    const int M = A->NZMAX_;
    double threshold = 0.0;
    bml_matrix_dense_t *D;
    bml_matrix_csr_t *B;
    REAL_T *A_dense;

    // From csr_bml to dense_array
    A_dense = bml_export_to_dense_csr(A, dense_row_major);

    // From dense_array to dense_bml
    D = bml_import_from_dense_dense(A->matrix_precision, dense_row_major,
                                    N, A_dense, threshold,
                                    A->distribution_mode);

    bml_free_memory(A_dense);

    // Calculate inverse in dense_bml
    TYPED_FUNC(bml_inverse_inplace_dense) (D);

    A_dense = bml_export_to_dense_dense(D, dense_row_major);
    bml_deallocate_dense(D);

    // From dense_array to csr_bml
    B = bml_import_from_dense_csr(A->matrix_precision,
                                  dense_row_major, N,
                                  A_dense,
                                  threshold, M, A->distribution_mode);
    bml_free_memory(A_dense);

    return B;
}
