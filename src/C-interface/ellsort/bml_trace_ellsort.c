#include "../bml_logger.h"
#include "../bml_trace.h"
#include "../bml_types.h"
#include "bml_trace_ellsort.h"
#include "bml_types_ellsort.h"

#include <stdlib.h>
#include <string.h>

/** Calculate the trace of a matrix.
 *
 *  \ingroup trace_group
 *
 *  \param A The matrix to calculate a trace for
 *  \return the trace of A
 */
double
bml_trace_ellsort(
    bml_matrix_ellsort_t * A)
{
    double trace = 0.0;

    switch (A->matrix_precision)
    {
        case single_real:
            trace = bml_trace_ellsort_single_real(A);
            break;
        case double_real:
            trace = bml_trace_ellsort_double_real(A);
            break;
        case single_complex:
            trace = bml_trace_ellsort_single_complex(A);
            break;
        case double_complex:
            trace = bml_trace_ellsort_double_complex(A);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
    return trace;
}

/** Calculate the trace of a matrix multiplication.
 * Both matrices must have the same size.
 *
 *  \ingroup trace_group
 *
 *  \param A The matrix A
 *  \param B The matrix B
 *  \return the trace of A*B
 */
double
bml_trace_mult_ellsort(
    bml_matrix_ellsort_t * A,
    bml_matrix_ellsort_t * B)
{
    double trace = 0.0;

    switch (A->matrix_precision)
    {
        case single_real:
            trace = bml_trace_mult_ellsort_single_real(A, B);
            break;
        case double_real:
            trace = bml_trace_mult_ellsort_double_real(A, B);
            break;
        case single_complex:
            trace = bml_trace_mult_ellsort_single_complex(A, B);
            break;
        case double_complex:
            trace = bml_trace_mult_ellsort_double_complex(A, B);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
    return trace;
}
