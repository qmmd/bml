#include "../macros.h"
#include "../typed.h"
#include "../C-interface/blas.h"
#include "../C-interface/bml_logger.h"

#include <complex.h>
#include <stdio.h>

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

    int m_val = *m;
    int n_val = *n;
    int k_val = *k;
    REAL_T alpha_val = *alpha;
    REAL_T beta_val = *beta;
 
    int N_rows_A;
    int N_rows_B;
    int N_cols_A;

    printf("alpha = %g, beta = %g\n",alpha_val,beta_val);
    
    if (*transa == 'N')
    {
        N_rows_A = m_val;
        N_cols_A = k_val;
    }
    else
    {
        N_rows_A = k_val;
        N_cols_A = m_val;
    }

    if (*transb == 'N')
    {
        N_rows_B = k_val;
    }
    else
    {
        N_rows_B = n_val;
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
    else if (m_val < 0)
    {
        info = 3;
    }
    else if (n_val < 0)
    {
        info = 4;
    }
    else if (k_val < 0)
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
    else if (*ldc < MAX(1, m_val))
    {
        info = 13;
    }

    if (info != 0)
    {
        /* Error. */
        LOG_ERROR("info = %d\n", info);
        return;
    }

    if ((m_val == 0 || n_val == 0) || ((alpha_val == 0 || k_val == 0) && beta_val == 1.0))
    {
        return;
    }

    if (alpha_val == 0)
    {
        if (beta_val == 0)
        {
	//#pragma omp target teams distribute parallel for simd collapse(2) schedule(static, 1)
	for (int j = 0; j < n_val; j++)
            {
                for (int i = 0; i < m_val; i++)
                {
                    c[COLMAJOR(i, j, m_val, n_val)] = 0;
                }
            }
        }
        else
        {
	//#pragma omp target teams distribute parallel for simd collapse(2) schedule(static, 1)
            for (int j = 0; j < n_val; j++)
            {
                for (int i = 0; i < m_val; i++)
                {
                    c[COLMAJOR(i, j, m_val, n_val)] =
                        beta_val * c[COLMAJOR(i, j, m_val, n_val)];
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

#pragma omp target teams distribute parallel for simd collapse(2) schedule(static, 1)
            for (int j = 0; j < n_val; j++)
            {
	    /*
                if (beta_val == 0)
                {
                    for (int i = 0; i < m_val; i++)
                    {
                        c[COLMAJOR(i, j, m_val, n_val)] = 0;
                    }
                }
                else if (beta_val != 1.0)
	    */
                {
                    for (int i = 0; i < m_val; i++)
                    {
                        c[COLMAJOR(i, j, m_val, n_val)] *= beta_val;
			for (int l = 0; l < k_val; l++) {
			c[COLMAJOR(i, j, m_val, n_val)] += alpha_val * a[COLMAJOR(i,l,m_val,k_val)] * b[COLMAJOR(l,j,k_val,n_val)];
			}
                    }
                }
	    }
	    
	    //#pragma omp target teams distribute parallel for simd collapse(2) schedule(static, 1)
	    /*            for (int j = 0; j < n_val; j++)
            {
	    for (int l = 0; l < k_val; l++)
                {
                    REAL_T temp = alpha_val * b[COLMAJOR(l, j, k_val, n_val)];
                    for (int i = 0; i < m_val; i++)
                    {
                        c[COLMAJOR(i, j, m_val, n_val)] +=
                            temp * a[COLMAJOR(i, l, m_val, k_val)];
                    }
                }
            }
	    */
        }
        else
        {
            /* C := alpha*A**T*B + beta*C
             */

	//#pragma omp target teams distribute parallel for simd collapse(2) schedule(static, 1)
            for (int j = 0; j < n_val; j++)
            {
                for (int i = 0; i < m_val; i++)
                {
                    REAL_T temp = 0;
                    for (int l = 0; l < k_val; l++)
                    {
                        temp +=
                            a[COLMAJOR(l, i, k_val, m_val)] *
                            b[COLMAJOR(l, j, k_val, n_val)];
                    }
		    /*
                    if (beta_val == 0)
                    {
                        c[COLMAJOR(i, j, m_val, n_val)] = alpha_val * temp;
                    }
                    else
		    */
                    {
                        c[COLMAJOR(i, j, m_val, n_val)] =
                            alpha_val * temp + beta_val * c[COLMAJOR(i, j, m_val, n_val)];
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

	//#pragma omp target teams distribute parallel for simd collapse(2) schedule(static, 1)
            for (int j = 0; j < n_val; ++j)
            {
	    /*
                if (beta_val == 0)
                {
                    for (int i = 0; i < m_val; i++)
                    {
                        c[COLMAJOR(i, j, m_val, n_val)] = 0;
                    }
                }
                else if (beta_val != 1.0)
                {
                    for (int i = 0; i < m_val; i++)
                    {
                        c[COLMAJOR(i, j, m_val, n_val)] *= beta_val;
                    }
                }

	    */
	    for (int l = 0; l < k_val; l++)
                {
                    REAL_T temp = alpha_val * b[COLMAJOR(j, l, n_val, k_val)];
                    for (int i = 0; i < m_val; i++)
                    {
                        c[COLMAJOR(i, j, m_val, n_val)] +=
                            temp * a[COLMAJOR(i, l, m_val, k_val)];
                    }
                }
            }
        }
        else
        {
            /* C := alpha*A**T*B**T + beta*C
             */

	//#pragma omp target teams distribute parallel for simd collapse(2) schedule(static, 1)
            for (int j = 0; j < n_val; j++)
            {
                for (int i = 0; i < m_val; i++)
                {
                    REAL_T temp = 0;
                    for (int l = 0; l < k_val; l++)
                    {
                        temp +=
                            a[COLMAJOR(l, i, k_val, m_val)] *
                            b[COLMAJOR(j, l, n_val, k_val)];
                    }

		    /*
                    if (beta_val == 0)
                    {
                        c[COLMAJOR(i, j, m_val, n_val)] = alpha_val * temp;
                    }
                    else
		    */
                    {
                        c[COLMAJOR(i, j, m_val, n_val)] =
                            alpha_val * temp + beta_val * c[COLMAJOR(i, j, m_val, n_val)];
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
