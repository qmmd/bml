#include "../macros.h"

#include <bml.h>
#include <math.h>
#include <stdlib.h>

//#define SINGLE_REAL
//#define DOUBLE_REAL
//#define SINGLE_COMPLEX
#define DOUBLE_COMPLEX

#ifdef SINGLE_REAL
#define REAL_T float
#define BML_T single_real
#define ABS fabsf
#elif defined(DOUBLE_REAL)
#define REAL_T double
#define BML_T double_real
#define ABS fabs
#elif defined(SINGLE_COMPLEX)
#define REAL_T complex float
#define BML_T single_complex
#define ABS cabsf
#elif defined(DOUBLE_COMPLEX)
#define REAL_T complex double
#define BML_T double_complex
#define ABS cabs
#endif

#define REL_TOL 1e-12
#define N 11
#define M 11

int main(int argc, char **argv)
{
    bml_matrix_t *A = NULL;
    REAL_T *A_dense = NULL;

    A_dense = calloc(N * N, sizeof(REAL_T));

    for (int i = 0; i < N * N; i++)
    {
        A_dense[i] = (REAL_T) (rand() / (REAL_T) RAND_MAX);
    }

    A = bml_import_from_dense(dense, BML_T, dense_row_major,
                              N, M, A_dense, 0.0, sequential);

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            REAL_T *Aij = bml_get(A, i, j);
            REAL_T expected = A_dense[ROWMAJOR(i, j, N, N)];
            double rel_diff = ABS((expected - *Aij) / expected);
            if (rel_diff > REL_TOL)
            {
                LOG_ERROR
                    ("matrices are not identical; expected[%d] = %e, B[%d] = %e\n",
                     i, expected, i, *Aij);
                return -1;
            }
        }
    }

    LOG_INFO("test passed\n");
}
