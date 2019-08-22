#include "../bml_logger.h"
#include "../bml_threshold.h"
#include "../bml_types.h"
#include "bml_threshold_dense.h"
#include "bml_types_dense.h"

#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/** Threshold a matrix.
 *
 *  \ingroup threshold_group
 *
 *  \param A The matrix to be thresholded
 *  \param threshold Threshold value
 *  \return The thresholded A
 */
bml_matrix_dense_t *
bml_threshold_new_dense(
    bml_matrix_dense_t * A,
    double threshold)
{
    switch (A->matrix_precision)
    {
        case single_real:
            return bml_threshold_new_dense_single_real(A, threshold);
            break;
        case double_real:
            return bml_threshold_new_dense_double_real(A, threshold);
            break;
        case single_complex:
            return bml_threshold_new_dense_single_complex(A, threshold);
            break;
        case double_complex:
            return bml_threshold_new_dense_double_complex(A, threshold);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
    return NULL;
}

/** Threshold a matrix in place.
 *
 *  \ingroup threshold_group
 *
 *  \param A The matrix to be thresholded
 *  \param threshold Threshold value
 *  \return The thresholded A
 */
void
bml_threshold_dense(
    bml_matrix_dense_t * A,
    double threshold)
{
    switch (A->matrix_precision)
    {
        case single_real:
            bml_threshold_dense_single_real(A, threshold);
            break;
        case double_real:
            bml_threshold_dense_double_real(A, threshold);
            break;
        case single_complex:
            bml_threshold_dense_single_complex(A, threshold);
            break;
        case double_complex:
            bml_threshold_dense_double_complex(A, threshold);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
}
