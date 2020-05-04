#include "bml_introspection.h"
#include "bml_logger.h"
#include "bml_normalize.h"
#include "dense/bml_normalize_dense.h"
#include "ellpack/bml_normalize_ellpack.h"
#include "ellsort/bml_normalize_ellsort.h"
#include "ellblock/bml_normalize_ellblock.h"
#include "csr/bml_normalize_csr.h"

#include <stdlib.h>

/** Normalize matrix given Gershgorin bounds.
 *
 ** \ingroup normalize_group_C
 *
 * \param A Matrix to scale
 * \param mineval Calculated min value
 * \param maxeval Calculated max value
 */
void
bml_normalize(
    bml_matrix_t * A,
    double mineval,
    double maxeval)
{
    switch (bml_get_type(A))
    {
        case dense:
            bml_normalize_dense(A, mineval, maxeval);
            break;
        case ellpack:
            bml_normalize_ellpack(A, mineval, maxeval);
            break;
        case ellsort:
            bml_normalize_ellsort(A, mineval, maxeval);
            break;
        case ellblock:
            bml_normalize_ellblock(A, mineval, maxeval);
            break;
        case csr:
            bml_normalize_csr(A, mineval, maxeval);
            break;
        default:
            LOG_ERROR("unknown matrix type\n");
            break;
    }
}

/** Calculate Gershgorin bounds.
 *
 * \ingroup normalize_group_C
 *
 * \param A Matrix to scale
 * returns mineval Calculated min value
 * returns maxeval Calculated max value
 */
void *
bml_gershgorin(
    bml_matrix_t * A)
{
    switch (bml_get_type(A))
    {
        case dense:
            return bml_gershgorin_dense(A);
            break;
        case ellpack:
            return bml_gershgorin_ellpack(A);
            break;
        case ellsort:
            return bml_gershgorin_ellsort(A);
            break;
        case ellblock:
            return bml_gershgorin_ellblock(A);
            break;
        case csr:
            return bml_gershgorin_csr(A);
            break;
        default:
            LOG_ERROR("unknown matrix type\n");
            break;
    }
    return NULL;
}

/** Calculate Gershgorin bounds for partial matrix.
 *
 * \ingroup normalize_group_C
 *
 * \param A Matrix to scale
 * \param nrows Number of rows used
 * returns mineval Calculated min value
 * returns maxeval Calculated max value
 */
void *
bml_gershgorin_partial(
    bml_matrix_t * A,
    int nrows)
{
    switch (bml_get_type(A))
    {
        case dense:
            return bml_gershgorin_partial_dense(A, nrows);
            break;
        case ellpack:
            return bml_gershgorin_partial_ellpack(A, nrows);
            break;
        case ellsort:
            return bml_gershgorin_partial_ellsort(A, nrows);
            break;
        case ellblock:
            LOG_ERROR("bml_gershgorin_partial_ellblock not implemented\n");
            break;
        case csr:
            return bml_gershgorin_partial_csr(A, nrows);
            break;
        default:
            LOG_ERROR("unknown matrix type\n");
            break;
    }
    return NULL;
}
