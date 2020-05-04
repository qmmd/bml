#include "bml_allocate_csr.h"
#include "../bml_allocate.h"
#include "bml_copy_csr.h"
#include "bml_diagonalize_csr.h"
#include "bml_export_csr.h"
#include "bml_import_csr.h"
#include "../bml_diagonalize.h"
#include "bml_types_csr.h"
#include "../bml_types.h"
#include "../bml_utilities.h"
#include "../dense/bml_allocate_dense.h"
#include "../dense/bml_import_dense.h"
#include "../dense/bml_diagonalize_dense.h"
#include "../dense/bml_types_dense.h"
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

/** Diagonalize matrix.
 *
 *  \ingroup diag_group
 *
 *  \param A The matrix A
 *  \param eigenvalues Eigenvalues of A
 *  \param eigenvectors Eigenvectors of A
 */
void TYPED_FUNC(
    bml_diagonalize_csr) (
    bml_matrix_csr_t * A,
    void *eigenvalues,
    bml_matrix_csr_t * eigenvectors)
{
    const int N = A->N_;
    const int M = A->NZMAX_;
    double threshold = 0.0;
    bml_matrix_dense_t *D;
    bml_matrix_dense_t *eigenvectors_bml_dense;
    REAL_T *A_dense;
    REAL_T *eigenvectors_dense;
    REAL_T *typed_eigenvalues = (REAL_T *) eigenvalues;
    bml_matrix_csr_t *myeigenvectors = (bml_matrix_csr_t *) eigenvectors;

    // Form csr_bml to dense_array
    A_dense = bml_export_to_dense_csr(A, dense_row_major);

    // From dense_array to dense_bml
    D = bml_import_from_dense_dense(A->matrix_precision, dense_row_major,
                                    N, A_dense, threshold,
                                    A->distribution_mode);
    free(A_dense);

    // Allocate eigenvectors matrix in dense_bml
    eigenvectors_bml_dense =
        bml_zero_matrix(dense, A->matrix_precision, N, N, sequential);

    // Diagonalize in dense_bml
    TYPED_FUNC(bml_diagonalize_dense) (D, typed_eigenvalues,
                                       eigenvectors_bml_dense);

    bml_deallocate_dense(D);

    // bml_print_bml_matrix(eigenvectors_bml_dense, 0, N, 0, N);

    // From dense_bml to dense_array
    eigenvectors_dense =
        bml_export_to_dense_dense(eigenvectors_bml_dense, dense_row_major);
    bml_deallocate_dense(eigenvectors_bml_dense);

    // From dense_array to csr_bml
    myeigenvectors = bml_import_from_dense_csr(A->matrix_precision,
                                               dense_row_major, N,
                                               eigenvectors_dense,
                                               threshold, M,
                                               A->distribution_mode);
    free(eigenvectors_dense);

    // This is done in order to pass the changes back to the upper level
    bml_copy_csr(myeigenvectors, eigenvectors);

    bml_deallocate_csr(myeigenvectors);

    return;
}
