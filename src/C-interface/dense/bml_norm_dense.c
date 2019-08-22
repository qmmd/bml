#include "../bml_logger.h"
#include "../bml_types.h"
#include "bml_norm_dense.h"
#include "bml_types_dense.h"

#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/** Calculate the sum of squares of the elements in matrix A.
 *
 *  \ingroup norm_group
 *
 *  \param A The matrix
 *  \return Sum of squares of all elements in A
 */
double
bml_sum_squares_dense(
    bml_matrix_dense_t * A)
{
    switch (A->matrix_precision)
    {
        case single_real:
            return bml_sum_squares_dense_single_real(A);
            break;
        case double_real:
            return bml_sum_squares_dense_double_real(A);
            break;
        case single_complex:
            return bml_sum_squares_dense_single_complex(A);
            break;
        case double_complex:
            return bml_sum_squares_dense_double_complex(A);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
    return 0;
}

/** Calculate the sum of squares of the core elements in submatrix A.
 *
 *  \ingroup norm_group
 *
 *  \param A The matrix
 *  \param core_pos Core rows in submatrix
 *  \param core_size Number of core rows
 *  \return Sum of squares of all elements in A
 */
double
bml_sum_squares_submatrix_dense(
    bml_matrix_dense_t * A,
    int core_size)
{
    switch (A->matrix_precision)
    {
        case single_real:
            return bml_sum_squares_submatrix_dense_single_real(A, core_size);
            break;
        case double_real:
            return bml_sum_squares_submatrix_dense_double_real(A, core_size);
            break;
        case single_complex:
            return bml_sum_squares_submatrix_dense_single_complex(A,
                                                                  core_size);
            break;
        case double_complex:
            return bml_sum_squares_submatrix_dense_double_complex(A,
                                                                  core_size);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
    return 0;
}

/** Calculate the sum of squares of elements in
 * \alpha A + \beta B.
 *
 *  \ingroup norm_group
 *
 *  \param A The matrix A
 *  \param B The matrix B
 *  \param alpha Multiplier for matrix A
 *  \param beta Multiplier for matrix B
 *  \param threshold Threshold
 *  \return The sum of squares of all elements of \alpha A + \beta B
 */
double
bml_sum_squares2_dense(
    bml_matrix_dense_t * A,
    bml_matrix_dense_t * B,
    double alpha,
    double beta,
    double threshold)
{
    switch (A->matrix_precision)
    {
        case single_real:
            return bml_sum_squares2_dense_single_real(A, B, alpha, beta,
                                                      threshold);
            break;
        case double_real:
            return bml_sum_squares2_dense_double_real(A, B, alpha, beta,
                                                      threshold);
            break;
        case single_complex:
            return bml_sum_squares2_dense_single_complex(A, B, alpha, beta,
                                                         threshold);
            break;
        case double_complex:
            return bml_sum_squares2_dense_double_complex(A, B, alpha, beta,
                                                         threshold);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
    return 0;
}

/** Calculate the Fobenius norm of matrix A.
 *
 *  \ingroup norm_group
 *
 *  \param A The matrix
 *  \return Frobenius norm of A
 */
double
bml_fnorm_dense(
    bml_matrix_dense_t * A)
{
    switch (A->matrix_precision)
    {
        case single_real:
            return bml_fnorm_dense_single_real(A);
            break;
        case double_real:
            return bml_fnorm_dense_double_real(A);
            break;
        case single_complex:
            return bml_fnorm_dense_single_complex(A);
            break;
        case double_complex:
            return bml_fnorm_dense_double_complex(A);
            break;
        default:
            LOG_ERROR("unknown precision");
            break;
    }
    return 0;
}

/** Calculate the Fobenius norm of 2 matrices.
 *
 *  \ingroup norm_group
 *
 *  \param A The matrix A
 *  \param B The matrix B
 *  \return Frobenius norm of A-B
 */
double
bml_fnorm2_dense(
    bml_matrix_dense_t * A,
    bml_matrix_dense_t * B)
{
    switch (A->matrix_precision)
    {
        case single_real:
            return bml_fnorm2_dense_single_real(A, B);
            break;
        case double_real:
            return bml_fnorm2_dense_double_real(A, B);
            break;
        case single_complex:
            return bml_fnorm2_dense_single_complex(A, B);
            break;
        case double_complex:
            return bml_fnorm2_dense_double_complex(A, B);
            break;
        default:
            LOG_ERROR("unknown precision");
            break;
    }
    return 0;
}
