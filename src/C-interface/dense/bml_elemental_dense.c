#include "../../macros.h"
#include "../bml_logger.h"
#include "bml_elemental_dense.h"
#include "bml_types_dense.h"

#include <complex.h>

/** Return a single matrix element.
 *
 * \param A The bml matrix
 * \param i The row index
 * \param j The column index
 * \return The matrix element
 */
float
bml_get_dense_single_real(
    const bml_matrix_dense_t * A,
    const int i,
    const int j)
{
    if (i < 0 || i >= A->N)
    {
        LOG_ERROR("row index out of bounds\n");
        return -1;
    }
    if (j < 0 || j >= A->N)
    {
        LOG_ERROR("column index out of bounds\n");
        return -1;
    }
    return ((float *) A->matrix)[LINEARINDEX(i, j, A->N, A->N)];
}

/** Return a single matrix element.
 *
 * \param A The bml matrix
 * \param i The row index
 * \param j The column index
 * \return The matrix element
 */
double
bml_get_dense_double_real(
    const bml_matrix_dense_t * A,
    const int i,
    const int j)
{
    if (i < 0 || i >= A->N)
    {
        LOG_ERROR("row index out of bounds\n");
        return -1;
    }
    if (j < 0 || j >= A->N)
    {
        LOG_ERROR("column index out of bounds\n");
        return -1;
    }
    return ((double *) A->matrix)[LINEARINDEX(i, j, A->N, A->N)];
}

/** Return a single matrix element.
 *
 * \param A The bml matrix
 * \param i The row index
 * \param j The column index
 * \return The matrix element
 */
float complex
bml_get_dense_single_complex(
    const bml_matrix_dense_t * A,
    const int i,
    const int j)
{
    if (i < 0 || i >= A->N)
    {
        LOG_ERROR("row index out of bounds\n");
        return -1;
    }
    if (j < 0 || j >= A->N)
    {
        LOG_ERROR("column index out of bounds\n");
        return -1;
    }
    return ((float complex *) A->matrix)[LINEARINDEX(i, j, A->N, A->N)];
}

/** Return a single matrix element.
 *
 * \param A The bml matrix
 * \param i The row index
 * \param j The column index
 * \return The matrix element
 */
double complex
bml_get_dense_double_complex(
    const bml_matrix_dense_t * A,
    const int i,
    const int j)
{
    if (i < 0 || i >= A->N)
    {
        LOG_ERROR("row index out of bounds\n");
        return -1;
    }
    if (j < 0 || j >= A->N)
    {
        LOG_ERROR("column index out of bounds\n");
        return -1;
    }
    return ((double complex *) A->matrix)[LINEARINDEX(i, j, A->N, A->N)];
}
