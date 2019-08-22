#include "../bml_utilities.h"
#include "../bml_logger.h"
#include "bml_types_ellsort.h"
#include "bml_utilities_ellsort.h"

void
bml_read_bml_matrix_ellsort(
    bml_matrix_ellsort_t * A,
    char *filename)
{
    switch (A->matrix_precision)
    {
        case single_real:
            bml_read_bml_matrix_ellsort_single_real(A, filename);
            break;
        case double_real:
            bml_read_bml_matrix_ellsort_double_real(A, filename);
            break;
        case single_complex:
            bml_read_bml_matrix_ellsort_single_complex(A, filename);
            break;
        case double_complex:
            bml_read_bml_matrix_ellsort_double_complex(A, filename);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
}

void
bml_write_bml_matrix_ellsort(
    bml_matrix_ellsort_t * A,
    char *filename)
{
    switch (A->matrix_precision)
    {
        case single_real:
            bml_write_bml_matrix_ellsort_single_real(A, filename);
            break;
        case double_real:
            bml_write_bml_matrix_ellsort_double_real(A, filename);
            break;
        case single_complex:
            bml_write_bml_matrix_ellsort_single_complex(A, filename);
            break;
        case double_complex:
            bml_write_bml_matrix_ellsort_double_complex(A, filename);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
}
