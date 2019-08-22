#include "../bml_utilities.h"
#include "../bml_logger.h"
#include "bml_types_dense.h"
#include "bml_utilities_dense.h"

void
bml_read_bml_matrix_dense(
    bml_matrix_dense_t * A,
    char *filename)
{
    switch (A->matrix_precision)
    {
        case single_real:
            bml_read_bml_matrix_dense_single_real(A, filename);
            break;
        case double_real:
            bml_read_bml_matrix_dense_double_real(A, filename);
            break;
        case single_complex:
            bml_read_bml_matrix_dense_single_complex(A, filename);
            break;
        case double_complex:
            bml_read_bml_matrix_dense_double_complex(A, filename);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
}

void
bml_write_bml_matrix_dense(
    bml_matrix_dense_t * A,
    char *filename)
{
    switch (A->matrix_precision)
    {
        case single_real:
            bml_write_bml_matrix_dense_single_real(A, filename);
            break;
        case double_real:
            bml_write_bml_matrix_dense_double_real(A, filename);
            break;
        case single_complex:
            bml_write_bml_matrix_dense_single_complex(A, filename);
            break;
        case double_complex:
            bml_write_bml_matrix_dense_double_complex(A, filename);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
}

void
bml_print_bml_matrix_dense(
    bml_matrix_dense_t * A,
    int i_l,
    int i_u,
    int j_l,
    int j_u)
{
    switch (A->matrix_precision)
    {
        case single_real:
            bml_print_bml_matrix_dense_single_real(A, i_l, i_u, j_l, j_u);
            break;
        case double_real:
            bml_print_bml_matrix_dense_double_real(A, i_l, i_u, j_l, j_u);
            break;
        case single_complex:
            bml_print_bml_matrix_dense_single_complex(A, i_l, i_u, j_l, j_u);
            break;
        case double_complex:
            bml_print_bml_matrix_dense_double_complex(A, i_l, i_u, j_l, j_u);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
}
