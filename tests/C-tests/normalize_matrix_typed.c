#include "bml.h"
#include "../typed.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>

#if defined(SINGLE_REAL) || defined(SINGLE_COMPLEX)
#define REL_TOL 1e-6
#else
#define REL_TOL 1e-12
#endif

int TYPED_FUNC(
    test_normalize) (
    const int N,
    const bml_matrix_type_t matrix_type,
    const bml_matrix_precision_t matrix_precision,
    const int M)
{
    bml_matrix_t *A = NULL;
    bml_matrix_t *B = NULL;
    double *A_gbnd = NULL;
    double *B_gbnd = NULL;
    REAL_T *A_dense = NULL;
    REAL_T *B_dense = NULL;

    REAL_T scale_factor = 2.5;
    double threshold = 0.0;

    A = bml_identity_matrix(matrix_type, matrix_precision, N, M, sequential);
    bml_scale_inplace(&scale_factor, A);
    A_gbnd = bml_gershgorin(A);

    A_dense = bml_export_to_dense(A, dense_row_major);
    A_dense[0] = scale_factor * scale_factor;
    B = bml_import_from_dense(matrix_type, matrix_precision, dense_row_major,
                              N, M, A_dense, threshold, sequential);
    B_gbnd = bml_gershgorin(B);

    bml_free_memory(A_dense);
    A_dense = bml_export_to_dense(A, dense_row_major);
    B_dense = bml_export_to_dense(B, dense_row_major);
    bml_print_dense_matrix(N, matrix_precision, dense_row_major, A_dense, 0,
                           N, 0, N);
    bml_print_dense_matrix(N, matrix_precision, dense_row_major, B_dense, 0,
                           N, 0, N);

    bml_normalize(B, B_gbnd[0], B_gbnd[1]);

    bml_free_memory(B_dense);
    B_dense = bml_export_to_dense(B, dense_row_major);
    bml_print_dense_matrix(N, matrix_precision, dense_row_major, B_dense, 0,
                           N, 0, N);

    if ((fabs(A_gbnd[1] - REAL_PART(scale_factor))) > REL_TOL
        || (A_gbnd[1] - A_gbnd[0]) > REL_TOL)
    {
        LOG_ERROR
            ("incorrect mineval or maxeval or maxminusmin; mineval = %e maxeval = %e maxminusmin = %e\n",
             A_gbnd[0], A_gbnd[1], A_gbnd[1] - A_gbnd[0]);
        return -1;
    }

    if ((fabs(B_gbnd[1] - REAL_PART(scale_factor * scale_factor))) > REL_TOL
        ||
        (fabs
         ((B_gbnd[1] - B_gbnd[0]) -
          REAL_PART(scale_factor * scale_factor - scale_factor))) > REL_TOL)
    {
        LOG_ERROR
            ("incorrect mineval or maxeval or maxminusmin; mineval = %e maxeval = %e maxminusmin = %e\n",
             B_gbnd[0], B_gbnd[1], B_gbnd[1] - B_gbnd[0]);
        return -1;
    }

    if (ABS(B_dense[0]) > REL_TOL)
    {
        LOG_ERROR
            ("normalize error, incorrect mineval or maxeval or maxminusmin; mineval = %e maxeval = %e maxminusmin = %e\n",
             B_gbnd[0], B_gbnd[1], B_gbnd[1] - B_gbnd[0]);
        return -1;
    }

    LOG_INFO("normalize matrix test passed\n");
    bml_free_memory(A_dense);
    bml_free_memory(B_dense);
    bml_free_memory(A_gbnd);
    bml_free_memory(B_gbnd);
    bml_deallocate(&A);
    bml_deallocate(&B);

    return 0;
}
