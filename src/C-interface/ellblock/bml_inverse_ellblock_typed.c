#include "bml_allocate_ellblock.h"
#include "../bml_allocate.h"
#include "bml_copy_ellblock.h"
#include "bml_export_ellblock.h"
#include "bml_inverse_ellblock.h"
#include "../bml_inverse.h"
#include "bml_types_ellblock.h"
#include "../bml_types.h"
#include "../bml_utilities.h"
#include "../dense/bml_allocate_dense.h"
#include "../dense/bml_import_dense.h"
#include "../dense/bml_inverse_dense.h"
#include "../dense/bml_types_dense.h"
#include "bml_import_ellblock.h"
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
bml_matrix_ellblock_t *TYPED_FUNC(
    bml_inverse_ellblock) (
    bml_matrix_ellblock_t * A)
{
    double threshold = 0.0;
    bml_matrix_dense_t *D;
    bml_matrix_ellblock_t *B;
    REAL_T *A_dense;

    // From ellblock_bml to dense_array
    A_dense = bml_export_to_dense_ellblock(A, dense_row_major);

    // From dense_array to dense_bml
    D = bml_import_from_dense_dense(A->matrix_precision, dense_row_major,
                                    A->N, A_dense, threshold,
                                    A->distribution_mode);

    bml_free_memory(A_dense);

    // Calculate inverse in dense_bml
    TYPED_FUNC(bml_inverse_inplace_dense) (D);

    A_dense = bml_export_to_dense_dense(D, dense_row_major);
    bml_deallocate_dense(D);

    // From dense_array to ellblock_bml
    B = bml_import_from_dense_ellblock(A->matrix_precision,
                                       dense_row_major, A->N,
                                       A_dense,
                                       threshold, A->M, A->distribution_mode);
    bml_free_memory(A_dense);

    return B;
}
