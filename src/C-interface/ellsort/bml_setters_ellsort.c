#include "../bml_introspection.h"
#include "../bml_logger.h"
#include "bml_setters_ellsort.h"
#include "bml_types_ellsort.h"


void
bml_set_element_new_ellsort(
    bml_matrix_ellsort_t * A,
    int i,
    int j,
    void *value)
{
    switch (bml_get_precision(A))
    {
        case single_real:
            bml_set_element_new_ellsort_single_real(A, i, j, value);
            break;
        case double_real:
            bml_set_element_new_ellsort_double_real(A, i, j, value);
            break;
        case single_complex:
            bml_set_element_new_ellsort_single_complex(A, i, j, value);
            break;
        case double_complex:
            bml_set_element_new_ellsort_double_complex(A, i, j, value);
            break;
        default:
            LOG_ERROR("unkonwn precision in bml_set_element_new_ellsort\n");
            break;
    }
}


void
bml_set_element_ellsort(
    bml_matrix_ellsort_t * A,
    int i,
    int j,
    void *value)
{
    switch (bml_get_precision(A))
    {
        case single_real:
            bml_set_element_ellsort_single_real(A, i, j, value);
            break;
        case double_real:
            bml_set_element_ellsort_double_real(A, i, j, value);
            break;
        case single_complex:
            bml_set_element_ellsort_single_complex(A, i, j, value);
            break;
        case double_complex:
            bml_set_element_ellsort_double_complex(A, i, j, value);
            break;
        default:
            LOG_ERROR("unkonwn precision in bml_set_element_ellsort\n");
            break;
    }
}

void
bml_set_row_ellsort(
    bml_matrix_ellsort_t * A,
    int i,
    void *row,
    double threshold)
{
    switch (A->matrix_precision)
    {
        case single_real:
            bml_set_row_ellsort_single_real(A, i, row, threshold);
            break;
        case double_real:
            bml_set_row_ellsort_double_real(A, i, row, threshold);
            break;
        case single_complex:
            bml_set_row_ellsort_single_complex(A, i, row, threshold);
            break;
        case double_complex:
            bml_set_row_ellsort_double_complex(A, i, row, threshold);
            break;
        default:
            LOG_ERROR("unkonwn precision\n");
            break;
    }
}

void
bml_set_diagonal_ellsort(
    bml_matrix_ellsort_t * A,
    void *diagonal,
    double threshold)
{
    switch (A->matrix_precision)
    {
        case single_real:
            bml_set_diagonal_ellsort_single_real(A, diagonal, threshold);
            break;
        case double_real:
            bml_set_diagonal_ellsort_double_real(A, diagonal, threshold);
            break;
        case single_complex:
            bml_set_diagonal_ellsort_single_complex(A, diagonal, threshold);
            break;
        case double_complex:
            bml_set_diagonal_ellsort_double_complex(A, diagonal, threshold);
            break;
        default:
            LOG_ERROR("unkonwn precision\n");
            break;
    }
}
