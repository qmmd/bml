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
    test_submatrix) (
    const int N,
    const bml_matrix_type_t matrix_type,
    const bml_matrix_precision_t matrix_precision,
    const int M)
{
    bml_matrix_t *A = NULL;
    bml_matrix_t *B = NULL;
    bml_matrix_t *C = NULL;
    bml_matrix_t *D = NULL;
    REAL_T *A_dense = NULL;
    REAL_T *B_dense = NULL;
    REAL_T *C_dense = NULL;
    REAL_T *D_dense = NULL;

    int nlist[1];
    int chlist[N];
    int vsize[2];
    int chlist_size, cpos_size;

    double threshold = 0.0;
    double threshold2 = 0.2;
    double threshold3 = 0.3;

    if (matrix_type == dense || matrix_type == ellblock || matrix_type == csr)
    {
        LOG_INFO("submatrix matrix test not available\n");
        return 0;
    }

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

    bml_threshold(A, threshold2);
    B = bml_threshold_new(A, threshold3);

    A_dense = bml_export_to_dense(A, dense_row_major);
    B_dense = bml_export_to_dense(B, dense_row_major);

    LOG_INFO("Matrix A:\n");
    bml_print_dense_matrix(N, matrix_precision, dense_row_major, A_dense, 0,
                           N, 0, N);
    LOG_INFO("Matrix B:\n");
    bml_print_dense_matrix(N, matrix_precision, dense_row_major, B_dense, 0,
                           N, 0, N);

    D = bml_zero_matrix(matrix_type, matrix_precision, N, M, sequential);

    for (int i = 0; i < N; i++)
    {
        // Get indices for submatrix
        nlist[0] = i;
        bml_matrix2submatrix_index(B, A, nlist, 1, chlist, vsize, 0);
        chlist_size = vsize[0];
        cpos_size = vsize[1];
        LOG_INFO("chlist size = %d cpos size %d\n", chlist_size, cpos_size);

        C = bml_zero_matrix(dense, matrix_precision, chlist_size,
                            chlist_size, sequential);

        // Get submatrix
        bml_matrix2submatrix(A, C, chlist, chlist_size);
        C_dense = bml_export_to_dense(C, dense_row_major);
        LOG_INFO("Submatrix C %d\n", i);
        bml_print_dense_matrix(chlist_size, matrix_precision, dense_row_major,
                               C_dense, 0, chlist_size, 0, chlist_size);
        // Collect into matrix
        bml_submatrix2matrix(C, D, chlist, chlist_size, cpos_size, threshold);
        D_dense = bml_export_to_dense(D, dense_row_major);
        LOG_INFO("Matrix D %d:\n", i);
        bml_print_dense_matrix(N, matrix_precision, dense_row_major, D_dense,
                               0, N, 0, N);

        bml_free_memory(C_dense);
        bml_free_memory(D_dense);
        bml_deallocate(&C);
    }

    D_dense = bml_export_to_dense(D, dense_row_major);
    LOG_INFO("Matrix D:\n");
    bml_print_dense_matrix(N, matrix_precision, dense_row_major, D_dense, 0,
                           N, 0, N);

    for (int i = 0; i < N * N; i++)
    {
        if (ABS(A_dense[i] - D_dense[i]) > REL_TOL)
        {
            LOG_ERROR("matrices are not identical A[%d] = %e D[%d] = %e\n", i,
                      A_dense[i], i, D_dense[i]);
            return -1;
        }
    }

    LOG_INFO("submatrix matrix test passed\n");
    bml_free_memory(A_dense);
    bml_free_memory(B_dense);
    bml_free_memory(D_dense);
    bml_deallocate(&A);
    bml_deallocate(&B);
    bml_deallocate(&D);

    return 0;
}
