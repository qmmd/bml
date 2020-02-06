#include "bml.h"
#include "../typed.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#if defined(SINGLE_REAL) || defined(SINGLE_COMPLEX)
#define REL_TOL 1e-6
#else
#define REL_TOL 1e-12
#endif

int TYPED_FUNC(
    test_adjungate_triangle) (
    const int N,
    const bml_matrix_type_t matrix_type,
    const bml_matrix_precision_t matrix_precision,
    const int M)
{
    bml_matrix_t * A = NULL;
    REAL_T *A_dense = NULL;

    if (matrix_type == dense || matrix_type == ellsort
        || matrix_type == ellblock)
    {
        LOG_INFO("adjungate triangle matrix test not available\n");
        return 0;
    }

    A_dense = bml_allocate_memory(sizeof(REAL_T) * N * N);
    for (int i = 0; i < N; i++)
    {
        int rsize = i * N;
        // Diagonal
        A_dense[rsize + i] = rand() / (double) RAND_MAX;
        if (i == 0)
        {
            A_dense[rsize + i + 1] = rand() / (double) RAND_MAX;
            A_dense[rsize + N + i] = rand() / (double) RAND_MAX;
            A_dense[rsize + i + 2] = rand() / (double) RAND_MAX;
        }
        else if (i < N - 1)
        {
            A_dense[rsize + i + 1] = rand() / (double) RAND_MAX;
            A_dense[rsize + N + i] = rand() / (double) RAND_MAX;
        }
        else
        {
            A_dense[rsize + i - 1] = rand() / (double) RAND_MAX;
            A_dense[rsize - N + i] = rand() / (double) RAND_MAX;
        }
    }
    A = bml_import_from_dense(matrix_type, matrix_precision, dense_row_major,
                              N, M, A_dense, 0, sequential);

    int *xadj = malloc(sizeof(int) * (N + 1));
    int *adjncy = malloc(sizeof(int) * (N * M));
    bml_adjacency(A, xadj, adjncy, 0);

    bml_print_dense_matrix(N, matrix_precision, dense_row_major, A_dense, 0,
                           N, 0, N);

/*
    for (int i = 0; i < N; i++)
    {
        printf("Row %d starts at %d\n", i, xadj[i]);
        for (int j = xadj[i]; j < xadj[i+1]; j++)
        {
            printf("%d ", adjncy[j]);
        }
        printf("\n");
    }
    printf("Total non-zeroes = %d\n", xadj[N]);
*/

    int idx = 0;
    for (int i = 0; i < N; i++)
    {
        if (xadj[i] != idx)
        {
            LOG_ERROR
                ("xadj values off, expected xadj[%d] = %d, actual xadj[%d] = %d\n",
                 i, idx, i, xadj[i]);
            return -1;
        }

        if (i == 0)
        {
            for (int j = 1; j < 3; j++)
            {
                if (adjncy[idx] != j)
                {
                    LOG_ERROR
                        ("adjncy values off, expected adjncy[%d] = %d, actual adjncy[%d] = %d\n",
                         idx, j, idx, adjncy[idx]);
                    return -1;
                }
                idx++;
            }
        }
        else if (i < N - 1)
        {
            for (int j = i - 1; j < i + 2; j += 2)
            {
                if (adjncy[idx] != j)
                {
                    LOG_ERROR
                        ("adjncy values off, expected adjncy[%d] = %d, actual adjncy[%d] = %d\n",
                         idx, j, idx, adjncy[idx]);
                    return -1;
                }
                idx++;
            }
        }
        else
        {
            if (adjncy[idx] != N - 2)
            {
                LOG_ERROR
                    ("adjncy values off, expected adjncy[%d] = %d, actual adjncy[%d] = %d\n",
                     idx, N - 2, idx, adjncy[idx]);
                return -1;
            }
            idx++;
        }
    }

    if (xadj[N] != idx)
    {
        LOG_ERROR
            ("xadj values off, expected xadj[%d] = %d, actual xadj[%d] = %d\n",
             N, idx, N, xadj[N]);
        return -1;
    }

    free(xadj);
    free(adjncy);
    bml_free_memory(A_dense);
    bml_deallocate(&A);


    LOG_INFO("adjacency matrix test passed\n");

    return 0;
}
