#include "bml.h"
#include "../typed.h"

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#if defined(SINGLE_REAL) || defined(SINGLE_COMPLEX)
#define REL_TOL 1e-6
#else
#define REL_TOL 1e-12
#endif

int TYPED_FUNC(
    test_trace) (
    const int N,
    const bml_matrix_type_t matrix_type,
    const bml_matrix_precision_t matrix_precision,
    const int M)
{
    bml_matrix_t * A = NULL;
    bml_matrix_t *B = NULL;
    bml_matrix_t *C = NULL;

    REAL_T *A_dense = NULL;
    REAL_T *B_dense = NULL;
    REAL_T *C_dense = NULL;

    double traceA, traceB, traceC;
    REAL_T scalar = 0.8;

    double rel_diff;

    LOG_DEBUG("rel. tolerance = %e\n", REL_TOL);

    A = bml_identity_matrix(matrix_type, matrix_precision, N, M, sequential);
    B = bml_scale_new(&scalar, A);
    C = bml_scale_new(&scalar, B);

    traceA = bml_trace(A);
    traceB = bml_trace(B);
    traceC = bml_trace(C);

    A_dense = bml_export_to_dense(A, dense_row_major);
    B_dense = bml_export_to_dense(B, dense_row_major);
    C_dense = bml_export_to_dense(C, dense_row_major);

    bml_print_dense_matrix(N, matrix_precision, dense_row_major, A_dense, 0,
                           N, 0, N);
    bml_print_dense_matrix(N, matrix_precision, dense_row_major, B_dense, 0,
                           N, 0, N);
    bml_print_dense_matrix(N, matrix_precision, dense_row_major, C_dense, 0,
                           N, 0, N);

#if defined(SINGLE_REAL) || defined(DOUBLE_REAL)
    printf("traceA = %e (%e), diff. traceA = %e\n", traceA, (double) N,
           traceA - (double) N);
    printf("traceB = %e (%e), diff. traceB = %e\n", traceB, scalar * N,
           traceB - scalar * N);
    printf("traceC = %e (%e), diff. traceC = %e\n", traceC,
           scalar * scalar * N, traceC - scalar * scalar * N);
#else
    printf("traceA = %e (%e), diff. traceA = %e\n", traceA, (double) N,
           traceA - (double) N);
    printf("traceB = %e (%e + %e i), diff. traceB = %e\n", traceB,
           REAL_PART(scalar * N), IMAGINARY_PART(scalar * N),
           traceB - REAL_PART(scalar * N));
    printf("traceC = %e (%e), diff. traceC = %e\n", traceC,
           REAL_PART(scalar * COMPLEX_CONJUGATE(scalar) * N),
           traceC - REAL_PART(scalar * COMPLEX_CONJUGATE(scalar) * N));
#endif

    if ((rel_diff = fabs(traceA - (double) N) / (double) N) > REL_TOL)
    {
        LOG_ERROR
            ("traces are not correct; traceA = %e and not %e, rel.diff = %e\n",
             traceA, N, rel_diff);
        return -1;
    }
    if ((rel_diff =
         fabs(traceB - (double) (scalar * N)) / (double) (scalar * N)) >
        REL_TOL)
    {
        LOG_ERROR
            ("traces are not correct; traceB = %e and not %e, rel.diff = %e\n",
             traceB, (double) (scalar * N), rel_diff);
        return -1;
    }
    if ((rel_diff =
         fabs(traceC -
              (double) (scalar * scalar * N)) / (double) (scalar * N * N)) >
        REL_TOL)
    {
        LOG_ERROR
            ("traces are not correct; traceC = %e and not %e, rel.diff = %e\n",
             traceC, (double) (scalar * scalar * N), rel_diff);
        return -1;
    }
    bml_free_memory(A_dense);
    bml_free_memory(B_dense);
    bml_free_memory(C_dense);
    bml_deallocate(&A);
    bml_deallocate(&B);
    bml_deallocate(&C);

    // Test when more values that just diagonal
    A_dense = bml_allocate_memory(sizeof(REAL_T) * N * N);
    for (int i = 0; i < N * N; i++)
    {
        A_dense[i] = rand() / (double) RAND_MAX;
    }
    for (int i = 0; i < N; i++)
    {
        A_dense[i * N + i] = (REAL_T) 1.0;
    }
    A = bml_import_from_dense(matrix_type, matrix_precision, dense_row_major,
                              N, M, A_dense, 0, sequential);
    bml_free_memory(A_dense);
    B = bml_scale_new(&scalar, A);
    C = bml_scale_new(&scalar, B);

    traceA = bml_trace(A);
    traceB = bml_trace(B);
    traceC = bml_trace(C);

    A_dense = bml_export_to_dense(A, dense_row_major);
    B_dense = bml_export_to_dense(B, dense_row_major);
    C_dense = bml_export_to_dense(C, dense_row_major);

    bml_print_dense_matrix(N, matrix_precision, dense_row_major, A_dense, 0,
                           N, 0, N);
    bml_print_dense_matrix(N, matrix_precision, dense_row_major, B_dense, 0,
                           N, 0, N);
    bml_print_dense_matrix(N, matrix_precision, dense_row_major, C_dense, 0,
                           N, 0, N);

#if defined(SINGLE_REAL) || defined(DOUBLE_REAL)
    printf("traceA = %e (%e), diff. traceA = %e\n", traceA, (double) N,
           traceA - (double) N);
    printf("traceB = %e (%e), diff. traceB = %e\n", traceB, scalar * N,
           traceB - scalar * N);
    printf("traceC = %e (%e), diff. traceC = %e\n", traceC,
           scalar * scalar * N, traceC - scalar * scalar * N);
#else
    printf("traceA = %e (%e), diff. traceA = %e\n", traceA, (double) N,
           traceA - (double) N);
    printf("traceB = %e (%e + %e i), diff. traceB = %e\n", traceB,
           REAL_PART(scalar * N), IMAGINARY_PART(scalar * N),
           traceB - REAL_PART(scalar * N));
    printf("traceC = %e (%e), diff. traceC = %e\n", traceC,
           REAL_PART(scalar * COMPLEX_CONJUGATE(scalar) * N),
           traceC - REAL_PART(scalar * COMPLEX_CONJUGATE(scalar) * N));
#endif

    if ((rel_diff = fabs(traceA - (double) N) / (double) N) > REL_TOL)
    {
        LOG_ERROR
            ("traces are not correct; traceA = %e and not %e, rel.diff = %e\n",
             traceA, rel_diff);
        return -1;
    }
    if ((rel_diff =
         fabs(traceB - (double) (scalar * N)) / (double) (scalar * N)) >
        REL_TOL)
    {
        LOG_ERROR
            ("traces are not correct; traceB = %e and not %e, rel.diff = %e\n",
             traceB, (double) (scalar * N), rel_diff);
        return -1;
    }
    if ((rel_diff =
         fabs(traceC -
              (double) (scalar * scalar * N)) / (double) (scalar * N * N)) >
        REL_TOL)
    {
        LOG_ERROR
            ("traces are not correct; traceC = %e and not %e, rel.diff = %e\n",
             traceC, (double) (scalar * scalar * N), rel_diff);
        return -1;
    }
    bml_free_memory(A_dense);
    bml_free_memory(B_dense);
    bml_free_memory(C_dense);
    bml_deallocate(&A);
    bml_deallocate(&B);
    bml_deallocate(&C);

    LOG_INFO("trace matrix test passed\n");

    return 0;
}
