#include "bml.h"
#include "../typed.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>


#if defined(SINGLE_REAL) || defined(SINGLE_COMPLEX)
#define REL_TOL 1.2e-5
#else
#define REL_TOL 1e-11
#endif

int TYPED_FUNC(
    test_diagonalize) (
    const int N,
    const bml_matrix_type_t matrix_type,
    const bml_matrix_precision_t matrix_precision,
    const int M)
{
    bml_matrix_t *A = NULL;
    bml_matrix_t *A_t = NULL;
    REAL_T *eigenvalues = NULL;
    bml_matrix_t *eigenvectors = NULL;
    bml_matrix_t *aux = NULL;
    bml_matrix_t *aux1 = NULL;
    bml_matrix_t *aux2 = NULL;
    float fnorm;

    if (matrix_type != dense && matrix_type != ellpack)
    {
        LOG_INFO("diagonalization for this matrix type is not implemented\n");
        return 0;
    }

    LOG_DEBUG("rel. tolerance = %e\n", REL_TOL);

    A = bml_random_matrix(matrix_type, matrix_precision, N, M, sequential);

    A_t = bml_transpose_new(A);
    bml_add(A, A_t, 0.5, 0.5, 0.0);
    bml_print_bml_matrix(A, 0, N, 0, N);

    switch (matrix_precision)
    {
        case single_real:
            eigenvalues = calloc(N, sizeof(float));
            break;
        case double_real:
            eigenvalues = calloc(N, sizeof(double));
            break;
        case single_complex:
            eigenvalues = calloc(N, sizeof(float complex));
            break;
        case double_complex:
            eigenvalues = calloc(N, sizeof(double complex));
            break;
        default:
            LOG_DEBUG("matrix_precision is not set");
            break;
    }

    eigenvectors = bml_zero_matrix(matrix_type, matrix_precision,
                                   N, M, sequential);

    aux1 = bml_zero_matrix(matrix_type, matrix_precision, N, M, sequential);
    aux2 = bml_zero_matrix(matrix_type, matrix_precision, N, M, sequential);

    bml_print_bml_matrix(eigenvectors, 0, N, 0, N);

    bml_diagonalize(A, eigenvalues, eigenvectors);

    printf("%s\n", "eigenvectors");
    bml_print_bml_matrix(eigenvectors, 0, N, 0, N);

    aux = bml_transpose_new(eigenvectors);

    bml_set_diagonal(aux1, eigenvalues, 0.0);

    bml_multiply(aux1, aux, aux2, 1.0, 0.0, 0.0);       // D*C^t

    bml_multiply(eigenvectors, aux2, aux, 1.0, 0.0, 0.0);       // C*(D*C^t)

    bml_add(aux, A, 1.0, -1.0, 0.0);

    fnorm = bml_fnorm(aux);

    if (fabsf(fnorm) > REL_TOL)
    {
        LOG_ERROR
            ("Error in matrix diagonalization; fnorm(CDC^t-A) = %e\n", fnorm);
        return -1;
    }

    bml_deallocate(&A);
    bml_deallocate(&aux);
    bml_deallocate(&aux1);
    bml_deallocate(&aux2);
    bml_deallocate(&A_t);
    bml_deallocate(&eigenvectors);
    free(eigenvalues);

    LOG_INFO("diagonalize matrix test passed\n");

    return 0;
}
