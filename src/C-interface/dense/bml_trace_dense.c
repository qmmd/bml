#include "../bml_logger.h"
#include "../bml_trace.h"
#include "../bml_types.h"
#include "bml_trace_dense.h"
#include "bml_types_dense.h"

#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/** Calculate the trace of a matrix.
 *
 *  \ingroup trace_group
 *
 *  \param A The matrix to calculate a trace for
 *  \return The trace of A
 */
double
bml_trace_dense(
    bml_matrix_dense_t * A)
{
    switch (A->matrix_precision)
    {
        case single_real:
            return bml_trace_dense_single_real(A);
            break;
        case double_real:
            return bml_trace_dense_double_real(A);
            break;
        case single_complex:
            return bml_trace_dense_single_complex(A);
            break;
        case double_complex:
            return bml_trace_dense_double_complex(A);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
    return 0;
}

/** Calculate the trace of a matrix multiplication.
 * Matrices must be of the same size.
 *
 *  \ingroup trace_group
 *
 *  \param A The matrix A
 *  \param B The matrix B
 *  \return The trace of A*B
 */
double
bml_trace_mult_dense(
    bml_matrix_dense_t * A,
    bml_matrix_dense_t * B)
{
    switch (A->matrix_precision)
    {
        case single_real:
            return bml_trace_mult_dense_single_real(A, B);
            break;
        case double_real:
            return bml_trace_mult_dense_double_real(A, B);
            break;
        case single_complex:
            return bml_trace_mult_dense_single_complex(A, B);
            break;
        case double_complex:
            return bml_trace_mult_dense_double_complex(A, B);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
    return 0;
}
