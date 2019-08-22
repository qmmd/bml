#include "../bml_introspection.h"
#include "../bml_logger.h"
#include "bml_setters_ellblock.h"
#include "bml_types_ellblock.h"


void
bml_set_element_new_ellblock(
    bml_matrix_ellblock_t * A,
    int i,
    int j,
    void *value)
{
    switch (bml_get_precision(A))
    {
        case single_real:
            bml_set_element_new_ellblock_single_real(A, i, j, value);
            break;
        case double_real:
            bml_set_element_new_ellblock_double_real(A, i, j, value);
            break;
        case single_complex:
            bml_set_element_new_ellblock_single_complex(A, i, j, value);
            break;
        case double_complex:
            bml_set_element_new_ellblock_double_complex(A, i, j, value);
            break;
        default:
            LOG_ERROR("unkonwn precision in bml_set_element_new_ellblock\n");
            break;
    }
}


void
bml_set_element_ellblock(
    bml_matrix_ellblock_t * A,
    int i,
    int j,
    void *value)
{
    switch (bml_get_precision(A))
    {
        case single_real:
            bml_set_element_ellblock_single_real(A, i, j, value);
            break;
        case double_real:
            bml_set_element_ellblock_double_real(A, i, j, value);
            break;
        case single_complex:
            bml_set_element_ellblock_single_complex(A, i, j, value);
            break;
        case double_complex:
            bml_set_element_ellblock_double_complex(A, i, j, value);
            break;
        default:
            LOG_ERROR("unkonwn precision in bml_set_element_ellblock\n");
            break;
    }
}

void
bml_set_row_ellblock(
    bml_matrix_ellblock_t * A,
    int i,
    void *row,
    double threshold)
{
    switch (A->matrix_precision)
    {
        case single_real:
            bml_set_row_ellblock_single_real(A, i, row, threshold);
            break;
        case double_real:
            bml_set_row_ellblock_double_real(A, i, row, threshold);
            break;
        case single_complex:
            bml_set_row_ellblock_single_complex(A, i, row, threshold);
            break;
        case double_complex:
            bml_set_row_ellblock_double_complex(A, i, row, threshold);
            break;
        default:
            LOG_ERROR("unkonwn precision\n");
            break;
    }
}

void
bml_set_diagonal_ellblock(
    bml_matrix_ellblock_t * A,
    void *diagonal,
    double threshold)
{
    switch (A->matrix_precision)
    {
        case single_real:
            bml_set_diagonal_ellblock_single_real(A, diagonal, threshold);
            break;
        case double_real:
            bml_set_diagonal_ellblock_double_real(A, diagonal, threshold);
            break;
        case single_complex:
            bml_set_diagonal_ellblock_single_complex(A, diagonal, threshold);
            break;
        case double_complex:
            bml_set_diagonal_ellblock_double_complex(A, diagonal, threshold);
            break;
        default:
            LOG_ERROR("unkonwn precision\n");
            break;
    }
}

void
bml_set_block_ellblock(
    bml_matrix_ellblock_t * A,
    int ib,
    int jb,
    void *values)
{
    switch (A->matrix_precision)
    {
        case single_real:
            bml_set_block_ellblock_single_real(A, ib, jb, values);
            break;
        case double_real:
            bml_set_block_ellblock_double_real(A, ib, jb, values);
            break;
        case single_complex:
            bml_set_block_ellblock_single_complex(A, ib, jb, values);
            break;
        case double_complex:
            bml_set_block_ellblock_double_complex(A, ib, jb, values);
            break;
        default:
            LOG_ERROR("unkonwn precision\n");
            break;
    }
}
