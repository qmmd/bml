#include "bml.h"
#include "../typed.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>

#if defined(SINGLE_REAL) || defined(SINGLE_COMPLEX)
#define REL_TOL 1e-3
#else
#define REL_TOL 1e-6
#endif

int TYPED_FUNC(
    test_norm) (
    const int N,
    const bml_matrix_type_t matrix_type,
    const bml_matrix_precision_t matrix_precision,
    const int M)
{
    bml_matrix_t *A = NULL;
    bml_matrix_t *B = NULL;
    REAL_T *A_dense = NULL;
    REAL_T *B_dense = NULL;

    const double alpha = 1.2;
    const double beta = 0.8;
    const double threshold = 0.0;

    double sum = 0.0;
    double sum2 = 0.0;
    double fnorm = 0.0;
    double sqrt_sum = 0.0;
    double sqrt_sum2 = 0.0;

    A = bml_random_matrix(matrix_type, matrix_precision, N, M, sequential);
    B = bml_random_matrix(matrix_type, matrix_precision, N, M, sequential);

    bml_print_bml_matrix(A, 0, N, 0, N);
    bml_print_bml_matrix(B, 0, N, 0, N);

    A_dense = bml_export_to_dense(A, dense_row_major);
    B_dense = bml_export_to_dense(B, dense_row_major);

    sum2 = bml_sum_squares2(A, B, alpha, beta, threshold);

    bml_add(A, B, alpha, beta, threshold);
    sum = bml_sum_squares(A);
    fnorm = bml_fnorm(A);

    //if (ABS(sum - sum2) > REL_TOL)
    if (fabs(sum - sum2) > REL_TOL)
    {
        LOG_ERROR
            ("incorrect sum of squares or sum of squares2; sum = %e sum2 = %e\n",
             sum, sum2);
        return -1;
    }

    sqrt_sum = sqrt(sum);
    sqrt_sum2 = sqrt(sum2);

    if ((fabs(sqrt_sum - sqrt_sum2) > REL_TOL) ||
        (fabs(sqrt_sum - fnorm) > REL_TOL) ||
        (fabs(sqrt_sum2 - fnorm) > REL_TOL))
    {
        LOG_ERROR
            ("incorrect sqrt(sum) or sqrt(sum2) of fnorm ; sqrt_sum = %e sqrt_sum2 = %e fnorm = %e\n",
             sqrt_sum, sqrt_sum2, fnorm);

        return -1;
    }

    LOG_INFO("norm matrix test passed\n");
    bml_free_memory(A_dense);
    bml_free_memory(B_dense);
    bml_deallocate(&A);
    bml_deallocate(&B);

    return 0;
}
