#include "bml_allocate_ellpack.h"
#include "../bml_allocate.h"
#include "bml_copy_ellpack.h"
#include "bml_diagonalize_ellpack.h"
#include "bml_export_ellpack.h"
#include "bml_import_ellpack.h"
#include "../bml_diagonalize.h"
#include "bml_types_ellpack.h"
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
    bml_diagonalize_ellpack) (
    const bml_matrix_ellpack_t * A,
    void *eigenvalues,
    bml_matrix_ellpack_t * eigenvectors)
{
    double threshold = 0.0;
    bml_matrix_dense_t *D;
    bml_matrix_dense_t *eigenvectors_bml_dense;
    REAL_T *A_dense;
    REAL_T *eigenvectors_dense;
    REAL_T *typed_eigenvalues = (REAL_T *) eigenvalues;
    bml_matrix_ellpack_t *myeigenvectors =
        (bml_matrix_ellpack_t *) eigenvectors;

    // Form ellpack_bml to dense_array
    A_dense = bml_export_to_dense_ellpack(A, dense_row_major);

    // From dense_array to dense_bml
    D = bml_import_from_dense_dense(A->matrix_precision, dense_row_major,
                                    A->N, A_dense, threshold,
                                    A->distribution_mode);
    free(A_dense);

    // Allocate eigenvectors matrix in dense_bml
    eigenvectors_bml_dense =
        bml_zero_matrix(dense, A->matrix_precision, A->N, A->N, sequential);

    // Diagonalize in dense_bml
    TYPED_FUNC(bml_diagonalize_dense) (D, typed_eigenvalues,
                                       eigenvectors_bml_dense);

    bml_deallocate_dense(D);

    // bml_print_bml_matrix(eigenvectors_bml_dense, 0, A->N, 0, A->N);

    // From dense_bml to dense_array
    eigenvectors_dense =
        bml_export_to_dense_dense(eigenvectors_bml_dense, dense_row_major);
    bml_deallocate_dense(eigenvectors_bml_dense);

    // From dense_array to ellpack_bml
    myeigenvectors = bml_import_from_dense_ellpack(A->matrix_precision,
                                                   dense_row_major, A->N,
                                                   eigenvectors_dense,
                                                   threshold, A->M,
                                                   A->distribution_mode);
    free(eigenvectors_dense);

    // This is done in order to pass the changes back to the upper level
    bml_copy_ellpack(myeigenvectors, eigenvectors);

    bml_deallocate_ellpack(myeigenvectors);

    return;
}
