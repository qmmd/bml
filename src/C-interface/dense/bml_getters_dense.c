#include "../bml_introspection.h"
#include "../bml_logger.h"
#include "bml_getters_dense.h"
#include "bml_types_dense.h"

void *
bml_get_dense(
    bml_matrix_dense_t * A,
    int i,
    int j)
{
    switch (bml_get_precision(A))
    {
        case single_real:
            return bml_get_dense_single_real(A, i, j);
            break;
        case double_real:
            return bml_get_dense_double_real(A, i, j);
            break;
        case single_complex:
            return bml_get_dense_single_complex(A, i, j);
            break;
        case double_complex:
            return bml_get_dense_double_complex(A, i, j);
            break;
        default:
            LOG_ERROR("unkonwn precision\n");
            break;
    }
    return NULL;
}

void *
bml_get_row_dense(
    bml_matrix_dense_t * A,
    int i)
{
    switch (bml_get_precision(A))
    {
        case single_real:
            return bml_get_row_dense_single_real(A, i);
            break;
        case double_real:
            return bml_get_row_dense_double_real(A, i);
            break;
        case single_complex:
            return bml_get_row_dense_single_complex(A, i);
            break;
        case double_complex:
            return bml_get_row_dense_double_complex(A, i);
            break;
        default:
            LOG_ERROR("unkonwn precision\n");
            break;
    }
    return NULL;
}

void *
bml_get_diagonal_dense(
    bml_matrix_dense_t * A)
{
    switch (bml_get_precision(A))
    {
        case single_real:
            return bml_get_diagonal_dense_single_real(A);
            break;
        case double_real:
            return bml_get_diagonal_dense_double_real(A);
            break;
        case single_complex:
            return bml_get_diagonal_dense_single_complex(A);
            break;
        case double_complex:
            return bml_get_diagonal_dense_double_complex(A);
            break;
        default:
            LOG_ERROR("unkonwn precision in bml_get_diagonal_dense\n");
            break;
    }
    return NULL;
}
