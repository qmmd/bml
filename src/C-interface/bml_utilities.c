#include "../macros.h"
#include "bml_export.h"
#include "bml_introspection.h"
#include "bml_logger.h"
#include "bml_types.h"
#include "bml_utilities.h"
#include "bml_allocate.h"
#include "dense/bml_utilities_dense.h"
#include "ellpack/bml_utilities_ellpack.h"
#include "ellsort/bml_utilities_ellsort.h"
#include "ellblock/bml_utilities_ellblock.h"
#include "csr/bml_utilities_csr.h"

#include <complex.h>
#include <stdio.h>

/** Print a bml vector.
 *
 * \param v The vector.
 * \param i_l The lower row index.
 * \param i_u The upper row index.
 */
void
bml_print_bml_vector(
    bml_vector_t * v,
    int i_l,
    int i_u)
{
    LOG_ERROR("[FIXME]\n");
}

/** Print a dense matrix.
 *
 * \param A The matrix.
 * \param i_l The lower row index.
 * \param i_u The upper row index.
 * \param j_l The lower column index.
 * \param j_u The upper column index.
 */
void
bml_print_bml_matrix(
    bml_matrix_t * A,
    int i_l,
    int i_u,
    int j_l,
    int j_u)
{
    switch (bml_get_type(A))
    {
        case dense:
            bml_print_bml_matrix_dense(A, i_l, i_u, j_l, j_u);
            break;
        case ellpack:
            switch (bml_get_precision(A))
            {
                case single_real:
                {
                    float *A_dense = bml_export_to_dense(A, dense_row_major);
                    bml_print_dense_matrix(bml_get_N(A), single_real,
                                           dense_row_major, A_dense, i_l, i_u,
                                           j_l, j_u);
                    bml_free_memory(A_dense);
                    break;
                }
                case double_real:
                {
                    double *A_dense = bml_export_to_dense(A, dense_row_major);
                    bml_print_dense_matrix(bml_get_N(A), double_real,
                                           dense_row_major, A_dense, i_l, i_u,
                                           j_l, j_u);
                    bml_free_memory(A_dense);
                    break;
                }
                case single_complex:
                {
                    float *A_dense = bml_export_to_dense(A, dense_row_major);
                    bml_print_dense_matrix(bml_get_N(A), single_complex,
                                           dense_row_major, A_dense, i_l, i_u,
                                           j_l, j_u);
                    bml_free_memory(A_dense);
                    break;
                }
                case double_complex:
                {
                    double *A_dense = bml_export_to_dense(A, dense_row_major);
                    bml_print_dense_matrix(bml_get_N(A), double_complex,
                                           dense_row_major, A_dense, i_l, i_u,
                                           j_l, j_u);
                    bml_free_memory(A_dense);
                    break;
                }
                default:
                    LOG_ERROR("unknown precision\n");
                    break;
            }
            break;
        case ellsort:
            switch (bml_get_precision(A))
            {
                case single_real:
                {
                    float *A_dense = bml_export_to_dense(A, dense_row_major);
                    bml_print_dense_matrix(bml_get_N(A), single_real,
                                           dense_row_major, A_dense, i_l, i_u,
                                           j_l, j_u);
                    bml_free_memory(A_dense);
                    break;
                }
                case double_real:
                {
                    double *A_dense = bml_export_to_dense(A, dense_row_major);
                    bml_print_dense_matrix(bml_get_N(A), double_real,
                                           dense_row_major, A_dense, i_l, i_u,
                                           j_l, j_u);
                    bml_free_memory(A_dense);
                    break;
                }
                case single_complex:
                {
                    float *A_dense = bml_export_to_dense(A, dense_row_major);
                    bml_print_dense_matrix(bml_get_N(A), single_complex,
                                           dense_row_major, A_dense, i_l, i_u,
                                           j_l, j_u);
                    bml_free_memory(A_dense);
                    break;
                }
                case double_complex:
                {
                    double *A_dense = bml_export_to_dense(A, dense_row_major);
                    bml_print_dense_matrix(bml_get_N(A), double_complex,
                                           dense_row_major, A_dense, i_l, i_u,
                                           j_l, j_u);
                    bml_free_memory(A_dense);
                    break;
                }
                default:
                    LOG_ERROR("unknown precision\n");
                    break;
            }
            break;
        case ellblock:
            switch (bml_get_precision(A))
            {
                case single_real:
                {
                    float *A_dense = bml_export_to_dense(A, dense_row_major);
                    bml_print_dense_matrix(bml_get_N(A), single_real,
                                           dense_row_major, A_dense, i_l, i_u,
                                           j_l, j_u);
                    bml_free_memory(A_dense);
                    break;
                }
                case double_real:
                {
                    double *A_dense = bml_export_to_dense(A, dense_row_major);
                    bml_print_dense_matrix(bml_get_N(A), double_real,
                                           dense_row_major, A_dense, i_l, i_u,
                                           j_l, j_u);
                    bml_free_memory(A_dense);
                    break;
                }
                case single_complex:
                {
                    float *A_dense = bml_export_to_dense(A, dense_row_major);
                    bml_print_dense_matrix(bml_get_N(A), single_complex,
                                           dense_row_major, A_dense, i_l, i_u,
                                           j_l, j_u);
                    bml_free_memory(A_dense);
                    break;
                }
                case double_complex:
                {
                    double *A_dense = bml_export_to_dense(A, dense_row_major);
                    bml_print_dense_matrix(bml_get_N(A), double_complex,
                                           dense_row_major, A_dense, i_l, i_u,
                                           j_l, j_u);
                    bml_free_memory(A_dense);
                    break;
                }
                default:
                    LOG_ERROR("unknown precision\n");
                    break;
            }
            break;
        case csr:
            switch (bml_get_precision(A))
            {
                case single_real:
                {
                    float *A_dense = bml_export_to_dense(A, dense_row_major);
                    bml_print_dense_matrix(bml_get_N(A), single_real,
                                           dense_row_major, A_dense, i_l, i_u,
                                           j_l, j_u);
                    bml_free_memory(A_dense);
                    break;
                }
                case double_real:
                {
                    double *A_dense = bml_export_to_dense(A, dense_row_major);
                    bml_print_dense_matrix(bml_get_N(A), double_real,
                                           dense_row_major, A_dense, i_l, i_u,
                                           j_l, j_u);
                    bml_free_memory(A_dense);
                    break;
                }
                case single_complex:
                {
                    float *A_dense = bml_export_to_dense(A, dense_row_major);
                    bml_print_dense_matrix(bml_get_N(A), single_complex,
                                           dense_row_major, A_dense, i_l, i_u,
                                           j_l, j_u);
                    bml_free_memory(A_dense);
                    break;
                }
                case double_complex:
                {
                    double *A_dense = bml_export_to_dense(A, dense_row_major);
                    bml_print_dense_matrix(bml_get_N(A), double_complex,
                                           dense_row_major, A_dense, i_l, i_u,
                                           j_l, j_u);
                    bml_free_memory(A_dense);
                    break;
                }
                default:
                    LOG_ERROR("unknown precision\n");
                    break;
            }
            break;
        default:
            LOG_ERROR("unknown type (%d)\n", bml_get_type(A));
            break;
    }
}

/** Print a dense matrix.
 *
 * \param N The number of rows/columns.
 * \param matrix_precision The real precision.
 * \param order The matrix element order.
 * \param A The matrix.
 * \param i_l The lower row index.
 * \param i_u The upper row index.
 * \param j_l The lower column index.
 * \param j_u The upper column index.
 */
void
bml_print_dense_matrix(
    int N,
    bml_matrix_precision_t matrix_precision,
    bml_dense_order_t order,
    void *A,
    int i_l,
    int i_u,
    int j_l,
    int j_u)
{
    if (N < 0)
    {
        LOG_ERROR("illegal value for N\n");
    }
    switch (matrix_precision)
    {
        case single_real:
        {
            float *A_typed = A;
            switch (order)
            {
                case dense_row_major:
                    for (int i = i_l; i < i_u; i++)
                    {
                        for (int j = j_l; j < j_u; j++)
                        {
                            printf(" % 1.3f", A_typed[ROWMAJOR(i, j, N, N)]);
                        }
                        printf("\n");
                    }
                    break;
                case dense_column_major:
                    for (int i = i_l; i < i_u; i++)
                    {
                        for (int j = j_l; j < j_u; j++)
                        {
                            printf(" % 1.3f", A_typed[COLMAJOR(i, j, N, N)]);
                        }
                        printf("\n");
                    }
                    break;
                default:
                    LOG_ERROR("logic error\n");
                    break;
            }
            break;
        }
        case double_real:
        {
            double *A_typed = A;
            switch (order)
            {
                case dense_row_major:
                    for (int i = i_l; i < i_u; i++)
                    {
                        for (int j = j_l; j < j_u; j++)
                        {
                            printf(" % 1.3f", A_typed[ROWMAJOR(i, j, N, N)]);
                        }
                        printf("\n");
                    }
                    break;
                case dense_column_major:
                    for (int i = i_l; i < i_u; i++)
                    {
                        for (int j = j_l; j < j_u; j++)
                        {
                            printf(" % 1.3f", A_typed[COLMAJOR(i, j, N, N)]);
                        }
                        printf("\n");
                    }
                    break;
                default:
                    LOG_ERROR("logic error\n");
                    break;
            }
            break;
        }
        case single_complex:
        {
            float complex *A_typed = A;
            switch (order)
            {
                case dense_row_major:
                    for (int i = i_l; i < i_u; i++)
                    {
                        for (int j = j_l; j < j_u; j++)
                        {
                            printf(" % 1.3f%+1.3fi",
                                   creal(A_typed[ROWMAJOR(i, j, N, N)]),
                                   cimag(A_typed[ROWMAJOR(i, j, N, N)]));
                        }
                        printf("\n");
                    }
                    break;
                case dense_column_major:
                    for (int i = i_l; i < i_u; i++)
                    {
                        for (int j = j_l; j < j_u; j++)
                        {
                            printf(" % 1.3f%+1.3fi",
                                   creal(A_typed[COLMAJOR(i, j, N, N)]),
                                   cimag(A_typed[COLMAJOR(i, j, N, N)]));
                        }
                        printf("\n");
                    }
                    break;
                default:
                    LOG_ERROR("logic error\n");
                    break;
            }
            break;
        }
        case double_complex:
        {
            double complex *A_typed = A;
            switch (order)
            {
                case dense_row_major:
                    for (int i = i_l; i < i_u; i++)
                    {
                        for (int j = j_l; j < j_u; j++)
                        {
                            printf(" % 1.3f%+1.3fi",
                                   creal(A_typed[ROWMAJOR(i, j, N, N)]),
                                   cimag(A_typed[ROWMAJOR(i, j, N, N)]));
                        }
                        printf("\n");
                    }
                    break;
                case dense_column_major:
                    for (int i = i_l; i < i_u; i++)
                    {
                        for (int j = j_l; j < j_u; j++)
                        {
                            printf(" % 1.3f%+1.3fi",
                                   creal(A_typed[COLMAJOR(i, j, N, N)]),
                                   cimag(A_typed[COLMAJOR(i, j, N, N)]));
                        }
                        printf("\n");
                    }
                    break;
                default:
                    LOG_ERROR("logic error\n");
                    break;
            }
            break;
        }
        default:
            LOG_ERROR("unknown matrix precision");
            break;
    }
}

/** Print a dense vector.
 *
 * \param N The number of rows/columns.
 * \param matrix_precision The real precision.
 * \param v The vector.
 * \param i_l The lower row index.
 * \param i_u The upper row index.
 */
void
bml_print_dense_vector(
    int N,
    bml_matrix_precision_t matrix_precision,
    void *v,
    int i_l,
    int i_u)
{
    switch (matrix_precision)
    {
        case single_real:
        {
            float *v_typed = v;
            for (int i = i_l; i < i_u; i++)
            {
                printf(" % 1.3f", v_typed[i]);
            }
            printf("\n");
            break;
        }
        case double_real:
        {
            double *v_typed = v;
            for (int i = i_l; i < i_u; i++)
            {
                printf(" % 1.3f", v_typed[i]);
            }
            printf("\n");
            break;
        }
        case single_complex:
        {
            float complex *v_typed = v;
            for (int i = i_l; i < i_u; i++)
            {
                printf(" % 1.3f%+1.3fi", creal(v_typed[i]),
                       cimag(v_typed[i]));
            }
            printf("\n");
            break;
        }
        case double_complex:
        {
            double complex *v_typed = v;
            for (int i = i_l; i < i_u; i++)
            {
                printf(" % 1.3f%+1.3fi", creal(v_typed[i]),
                       cimag(v_typed[i]));
            }
            printf("\n");
            break;
        }
        default:
            LOG_ERROR("unknown matrix precision");
            break;
    }
}

/** Read a bml matrix from a Matrix Market file.
 *
 * \param A The matrix
 * \param filename The file containing matrix
 */
void
bml_read_bml_matrix(
    bml_matrix_t * A,
    char *filename)
{
    switch (bml_get_type(A))
    {
        case dense:
            bml_read_bml_matrix_dense(A, filename);
            break;
        case ellpack:
            bml_read_bml_matrix_ellpack(A, filename);
            break;
        case ellsort:
            bml_read_bml_matrix_ellsort(A, filename);
            break;
        case ellblock:
            bml_read_bml_matrix_ellblock(A, filename);
            break;
        case csr:
            bml_read_bml_matrix_csr(A, filename);
            break;
        default:
            LOG_ERROR("unknown type (%d)\n", bml_get_type(A));
            break;
    }
}

/** Write a bml matrix to a Matrix Market file.
 *
 * \param A The matrix
 * \param filename The file containing matrix
 */
void
bml_write_bml_matrix(
    bml_matrix_t * A,
    char *filename)
{
    switch (bml_get_type(A))
    {
        case dense:
            bml_write_bml_matrix_dense(A, filename);
            break;
        case ellpack:
            bml_write_bml_matrix_ellpack(A, filename);
            break;
        case ellsort:
            bml_write_bml_matrix_ellsort(A, filename);
            break;
        case ellblock:
            bml_write_bml_matrix_ellblock(A, filename);
            break;
        case csr:
            bml_write_bml_matrix_csr(A, filename);
            break;
        default:
            LOG_ERROR("unknown type (%d)\n", bml_get_type(A));
            break;
    }
}

// compute largest int "i" such that i*i <= x
int
bml_sqrtint(
    const int x)
{
    int i = 1;
    int result = 1;
    while (result < x)
    {
        i++;
        result = i * i;
    }
    return i;
}
