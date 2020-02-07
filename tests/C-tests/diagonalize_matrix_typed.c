#include "bml.h"
#include "../typed.h"
#include "../C-interface/dense/bml_getters_dense.h"

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
    bml_matrix_t *id = NULL;
    float fnorm;

    LOG_DEBUG("rel. tolerance = %e\n", REL_TOL);

    A = bml_random_matrix(matrix_type, matrix_precision, N, M, sequential);

    A_t = bml_transpose_new(A);
    bml_add(A, A_t, 0.5, 0.5, 0.0);
    bml_print_bml_matrix(A, 0, N, 0, N);

    switch (matrix_precision)
    {
        case single_real:
#ifdef INTEL_OPT
            eigenvalues = bml_allocate_memory(N * sizeof(float));

#pragma omp parallel for simd
#pragma vector aligned
            for (int i = 0; i < N; i++)
            {
                __assume_aligned(eigenvalues, 64);
                eigenvalues[i] = 0.0;
            }
#else
            eigenvalues = calloc(N, sizeof(float));
#endif
            break;
        case double_real:
#ifdef INTEL_OPT
            eigenvalues = bml_allocate_memory(N * sizeof(double));
#pragma omp parallel for simd
#pragma vector aligned
            for (int i = 0; i < N; i++)
            {
                __assume_aligned(eigenvalues, 64);
                eigenvalues[i] = 0.0;
            }
#else
            eigenvalues = calloc(N, sizeof(double));
#endif
            break;
        case single_complex:
#ifdef INTEL_OPT
            eigenvalues = bml_allocate_memory(N * sizeof(float complex));
#pragma omp parallel for simd
#pragma vector aligned
            for (int i = 0; i < N; i++)
            {
                __assume_aligned(eigenvalues, 64);
                eigenvalues[i] = 0.0;
            }
#else
            eigenvalues = calloc(N, sizeof(float complex));
#endif
            break;
        case double_complex:
#ifdef INTEL_OPT
            eigenvalues = bml_allocate_memory(N * sizeof(double complex));
#pragma omp parallel for simd
#pragma vector aligned
            for (int i = 0; i < N; i++)
            {
                __assume_aligned(eigenvalues, 64);
                eigenvalues[i] = 0.0;
            }
#else
            eigenvalues = calloc(N, sizeof(double complex));
#endif
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

    printf("%s\n", "eigenvalues");
    for (int i = 0; i < N; i++)
        printf("val = %e  i%e\n", REAL_PART(eigenvalues[i]),
               IMAGINARY_PART(eigenvalues[i]));

    aux = bml_transpose_new(eigenvectors);
    bml_multiply(aux, eigenvectors, aux2, 1.0, 0.0, 0.0);       // C^t*C
    printf("%s\n", "check eigenvectors norms");
    for (int i = 0; i < N; i++)
    {
        REAL_T *val = bml_get(aux2, i, i);
        if (ABS(*val - (REAL_T) 1.0) > REL_TOL)
        {
            printf("i = %d, val = %e  i%e\n", i, REAL_PART(*val),
                   IMAGINARY_PART(*val));
            LOG_ERROR
                ("Error in matrix diagonalization; eigenvector not normalized\n");
        }
    }

    id = bml_identity_matrix(matrix_type, matrix_precision, N, M, sequential);
    bml_add(aux2, id, 1.0, -1.0, 0.0);
    fnorm = bml_fnorm(aux2);
    if (fabsf(fnorm) > N * REL_TOL)
    {
        LOG_ERROR
            ("Error in matrix diagonalization; fnorm(C^txC^t-Id) = %e\n",
             fnorm);
        return -1;
    }

    bml_set_diagonal(aux1, eigenvalues, 0.0);

    bml_multiply(aux1, aux, aux2, 1.0, 0.0, 0.0);       // D*C^t

    bml_multiply(eigenvectors, aux2, aux, 1.0, 0.0, 0.0);       // C*(D*C^t)

    bml_add(aux, A, 1.0, -1.0, 0.0);

    fnorm = bml_fnorm(aux);

    if (fabsf(fnorm) > N * REL_TOL)
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
    bml_deallocate(&id);
    bml_free_memory(eigenvalues);

    LOG_INFO("diagonalize matrix test passed\n");

    return 0;
}
