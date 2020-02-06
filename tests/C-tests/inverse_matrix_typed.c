#include "bml.h"
#include "../typed.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>


#if defined(SINGLE_REAL) || defined(SINGLE_COMPLEX)
#define REL_TOL 5e-4
#else
#define REL_TOL 1e-11
#endif

int TYPED_FUNC(
    test_inverse) (
    const int N,
    const bml_matrix_type_t matrix_type,
    const bml_matrix_precision_t matrix_precision,
    const int M)
{
    bml_matrix_t * A = NULL;
    bml_matrix_t *A_inverse = NULL;
    bml_matrix_t *aux = NULL;
    float ssum;

    LOG_DEBUG("rel. tolerance = %e\n", REL_TOL);

    A = bml_random_matrix(matrix_type, matrix_precision, N, M, sequential);

    A_inverse = bml_inverse(A);

    LOG_INFO("A:\n");
    bml_print_bml_matrix(A, 0, N, 0, N);
    LOG_INFO("A_inverse:\n");
    bml_print_bml_matrix(A_inverse, 0, N, 0, N);

    aux = bml_zero_matrix(matrix_type, matrix_precision, N, M, sequential);

    bml_multiply(A, A_inverse, aux, 1.0, 0.0, 0.0);     // A*A_inverse
    LOG_INFO("A*A^{-1}:\n");
    bml_print_bml_matrix(aux, 0, N, 0, N);

    ssum = bml_sum_squares(aux) - (float) N;

    if (fabsf(ssum) > REL_TOL)
    {
        LOG_ERROR("Error in matrix inverse; ssum(A*A_inverse) = %e\n", ssum);
        return -1;
    }

    bml_deallocate(&A);
    bml_deallocate(&A_inverse);
    bml_deallocate(&aux);

    LOG_INFO("inverse matrix test passed\n");

    return 0;
}
