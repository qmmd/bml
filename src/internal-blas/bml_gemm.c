#include "../macros.h"
#include "../typed.h"
#include "../C-interface/blas.h"
#include "../C-interface/bml_logger.h"

#include <complex.h>

void TYPED_FUNC(
    bml_gemm_internal) (
    const char *transa,
    const char *transb,
    const int *m,
    const int *n,
    const int *k,
    const REAL_T * alpha,
    const REAL_T * a,
    const int *lda,
    const REAL_T * b,
    const int *ldb,
    const REAL_T * beta,
    REAL_T * c,
    const int *ldc)
{
    /* Reference implementation from
     * http://www.netlib.org/lapack/explore-html/d1/d54/group__double__blas__level3_gaeda3cbd99c8fb834a60a6412878226e1.html
     */

    int N_rows_A;
    int N_rows_B;
    int N_cols_A;

    if (*transa == 'N')
    {
        N_rows_A = *m;
        N_cols_A = *k;
    }
    else
    {
        N_rows_A = *k;
        N_cols_A = *m;
    }

    if (*transb == 'N')
    {
        N_rows_B = *k;
    }
    else
    {
        N_rows_B = *n;
    }

    int info = 0;

    if (*transa != 'N' && *transa != 'C' && *transa != 'T')
    {
        info = 1;
    }
    else if (*transb != 'N' && *transb != 'C' && *transb != 'T')
    {
        info = 2;
    }
    else if (*m < 0)
    {
        info = 3;
    }
    else if (*n < 0)
    {
        info = 4;
    }
    else if (*k < 0)
    {
        info = 5;
    }
    else if (*lda < MAX(1, N_rows_A))
    {
        info = 8;
    }
    else if (*ldb < MAX(1, N_rows_B))
    {
        info = 10;
    }
    else if (*ldc < MAX(1, *m))
    {
        info = 13;
    }

    if (info != 0)
    {
        /* Error. */
        LOG_ERROR("info = %d\n", info);
        return;
    }

    if ((*m == 0 || *n == 0) || ((*alpha == 0 || *k == 0) && *beta == 1.0))
    {
        return;
    }

    if (*alpha == 0)
    {
        if (*beta == 0)
        {
            for (int j = 0; j < *n; j++)
            {
                for (int i = 0; i < *m; i++)
                {
                    c[COLMAJOR(i, j, *m, *n)] = 0;
                }
            }
        }
        else
        {
            for (int j = 0; j < *n; j++)
            {
                for (int i = 0; i < *m; i++)
                {
                    c[COLMAJOR(i, j, *m, *n)] =
                        *beta * c[COLMAJOR(i, j, *m, *n)];
                }
            }
        }
        return;
    }

    if (*transb == 'N')
    {
        if (*transa == 'N')
        {
            /* C := alpha*A*B + beta*C
             */

            for (int j = 0; j < *n; j++)
            {
                if (*beta == 0)
                {
                    for (int i = 0; i < *m; i++)
                    {
                        c[COLMAJOR(i, j, *m, *n)] = 0;
                    }
                }
                else if (*beta != 1.0)
                {
                    for (int i = 0; i < *m; i++)
                    {
                        c[COLMAJOR(i, j, *m, *n)] *= *beta;
                    }
                }

                for (int l = 0; l < *k; l++)
                {
                    REAL_T temp = *alpha * b[COLMAJOR(l, j, *k, *n)];
                    for (int i = 0; i < *m; i++)
                    {
                        c[COLMAJOR(i, j, *m, *n)] +=
                            temp * a[COLMAJOR(i, l, *m, *k)];
                    }
                }
            }
        }
        else
        {
            /* C := alpha*A**T*B + beta*C
             */

            for (int j = 0; j < *n; j++)
            {
                for (int i = 0; i < *m; i++)
                {
                    REAL_T temp = 0;
                    for (int l = 0; l < *k; l++)
                    {
                        temp +=
                            a[COLMAJOR(l, i, *k, *m)] *
                            b[COLMAJOR(l, j, *k, *n)];
                    }
                    if (*beta == 0)
                    {
                        c[COLMAJOR(i, j, *m, *n)] = *alpha * temp;
                    }
                    else
                    {
                        c[COLMAJOR(i, j, *m, *n)] =
                            *alpha * temp + *beta * c[COLMAJOR(i, j, *m, *n)];
                    }
                }
            }
        }
    }
    else
    {
        if (*transa == 'N')
        {
            /* C := alpha*A*B**T + beta*C
             */

            for (int j = 0; j < *n; ++j)
            {
                if (*beta == 0)
                {
                    for (int i = 0; i < *m; i++)
                    {
                        c[COLMAJOR(i, j, *m, *n)] = 0;
                    }
                }
                else if (*beta != 1.0)
                {
                    for (int i = 0; i < *m; i++)
                    {
                        c[COLMAJOR(i, j, *m, *n)] *= *beta;
                    }
                }

                for (int l = 0; l < *k; l++)
                {
                    REAL_T temp = *alpha * b[COLMAJOR(j, l, *n, *k)];
                    for (int i = 0; i < *m; i++)
                    {
                        c[COLMAJOR(i, j, *m, *n)] +=
                            temp * a[COLMAJOR(i, l, *m, *k)];
                    }
                }
            }
        }
        else
        {
            /* C := alpha*A**T*B**T + beta*C
             */

            for (int j = 0; j < *n; j++)
            {
                for (int i = 0; i < *m; i++)
                {
                    REAL_T temp = 0;
                    for (int l = 0; l < *k; l++)
                    {
                        temp +=
                            a[COLMAJOR(l, i, *k, *m)] *
                            b[COLMAJOR(j, l, *n, *k)];
                    }

                    if (*beta == 0)
                    {
                        c[COLMAJOR(i, j, *m, *n)] = *alpha * temp;
                    }
                    else
                    {
                        c[COLMAJOR(i, j, *m, *n)] =
                            *alpha * temp + *beta * c[COLMAJOR(i, j, *m, *n)];
                    }
                }
            }
        }
    }
}

void TYPED_FUNC(
    bml_gemm) (
    const char *transa,
    const char *transb,
    const int *m,
    const int *n,
    const int *k,
    const REAL_T * alpha,
    const REAL_T * a,
    const int *lda,
    const REAL_T * b,
    const int *ldb,
    const REAL_T * beta,
    REAL_T * c,
    const int *ldc)
{
#ifdef BML_INTERNAL_GEMM
    TYPED_FUNC(bml_gemm_internal) (transa, transb, m, n, k, alpha, a,
                                   lda, b, ldb, beta, c, ldc);
#else

#ifdef NOBLAS
    LOG_ERROR("No BLAS library");
#else
    C_BLAS(GEMM) (transa, transb, m, n, k, alpha, a,
                  lda, b, ldb, beta, c, ldc);
#endif
#endif
}
