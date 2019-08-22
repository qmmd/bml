#include "bml_introspection.h"
#include "bml_logger.h"
#include "bml_setters.h"
#include "dense/bml_setters_dense.h"
#include "ellpack/bml_setters_ellpack.h"
#include "ellsort/bml_setters_ellsort.h"
#include "ellblock/bml_setters_ellblock.h"

void
bml_set_element_new(
    bml_matrix_t * A,
    int i,
    int j,
    void *value)
{
    switch (bml_get_type(A))
    {
        case dense:
            bml_set_element_dense(A, i, j, value);
            break;
        case ellpack:
            bml_set_element_new_ellpack(A, i, j, value);
            break;
        case ellsort:
            bml_set_element_new_ellsort(A, i, j, value);
            break;
        case ellblock:
            bml_set_element_new_ellblock(A, i, j, value);
            break;
        default:
            LOG_ERROR("unknown matrix type in bml_set_new\n");
            break;
    }
}


void
bml_set_element(
    bml_matrix_t * A,
    int i,
    int j,
    void *value)
{
    switch (bml_get_type(A))
    {
        case dense:
            bml_set_element_dense(A, i, j, value);
            break;
        case ellpack:
            bml_set_element_ellpack(A, i, j, value);
            break;
        case ellsort:
            bml_set_element_ellsort(A, i, j, value);
            break;
        case ellblock:
            bml_set_element_ellblock(A, i, j, value);
            break;
        default:
            LOG_ERROR("unknown matrix type in bml_set\n");
            break;
    }
}

void
bml_set_row(
    bml_matrix_t * A,
    int i,
    void *row,
    double threshold)
{
    switch (bml_get_type(A))
    {
        case dense:
            bml_set_row_dense(A, i, row);
            break;
        case ellpack:
            bml_set_row_ellpack(A, i, row, threshold);
            break;
        case ellsort:
            bml_set_row_ellsort(A, i, row, threshold);
            break;
        case ellblock:
            bml_set_row_ellblock(A, i, row, threshold);
            break;
        default:
            LOG_ERROR("unknown matrix type in bml_set_row\n");
            break;
    }
}

void
bml_set_diagonal(
    bml_matrix_t * A,
    void *diagonal,
    double threshold)
{
    switch (bml_get_type(A))
    {
        case dense:
            bml_set_diagonal_dense(A, diagonal);
            break;
        case ellpack:
            bml_set_diagonal_ellpack(A, diagonal, threshold);
            break;
        case ellsort:
            bml_set_diagonal_ellsort(A, diagonal, threshold);
            break;
        case ellblock:
            bml_set_diagonal_ellblock(A, diagonal, threshold);
            break;
        default:
            LOG_ERROR("unknown matrix type in bml_set_diagonal\n");
            break;
    }
}
