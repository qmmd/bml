#include "bml_add.h"
#include "bml_introspection.h"
#include "bml_logger.h"
#include "dense/bml_add_dense.h"
#include "ellpack/bml_add_ellpack.h"
#include "ellsort/bml_add_ellsort.h"
#include "ellblock/bml_add_ellblock.h"

#include <stdlib.h>

/** Matrix addition.
 *
 * \f$ A \leftarrow \alpha A + \beta B \f$
 *
 * \ingroup add_group_C
 *
 * \param A Matrix A
 * \param B Matrix B
 * \param alpha Scalar factor multiplied by A
 * \param beta Scalar factor multiplied by B
 * \param threshold Threshold for matrix addition
 */
void
bml_add(
    bml_matrix_t * A,
    bml_matrix_t * B,
    double alpha,
    double beta,
    double threshold)
{
    switch (bml_get_type(A))
    {
        case dense:
            bml_add_dense(A, B, alpha, beta);
            break;
        case ellpack:
            bml_add_ellpack(A, B, alpha, beta, threshold);
            break;
        case ellsort:
            bml_add_ellsort(A, B, alpha, beta, threshold);
            break;
        case ellblock:
            bml_add_ellblock(A, B, alpha, beta, threshold);
            break;
        default:
            LOG_ERROR("unknown matrix type\n");
            break;
    }
}

/** Matrix addition with calculation of TrNorm.
 *
 * \f$ A \leftarrow \alpha A + \beta B \f$
 *
 * \ingroup add_group_C
 *
 * \param A Matrix A
 * \param B Matrix B
 * \param alpha Scalar factor multiplied by A
 * \param beta Scalar factor multiplied by B
 * \param threshold Threshold for matrix addition
 *
 */
double
bml_add_norm(
    bml_matrix_t * A,
    bml_matrix_t * B,
    double alpha,
    double beta,
    double threshold)
{

    switch (bml_get_type(A))
    {
        case dense:
            return bml_add_norm_dense(A, B, alpha, beta);
            break;
        case ellpack:
            return bml_add_norm_ellpack(A, B, alpha, beta, threshold);
            break;
        case ellblock:
            return bml_add_norm_ellblock(A, B, alpha, beta, threshold);
            break;
        case ellsort:
            return bml_add_norm_ellsort(A, B, alpha, beta, threshold);
            break;
        default:
            LOG_ERROR("unknown matrix type\n");
            break;
    }
    return 0;
}

/** Matrix addition.
 *
 * \f$ A \leftarrow A + \beta \mathrm{Id} \f$
 *
 * \ingroup add_group_C
 *
 * \param A Matrix A
 * \param beta Scalar factor multiplied by I
 * \param threshold Threshold for matrix addition
 */
void
bml_add_identity(
    bml_matrix_t * A,
    double beta,
    double threshold)
{
    switch (bml_get_type(A))
    {
        case dense:
            bml_add_identity_dense(A, beta);
            break;
        case ellpack:
            bml_add_identity_ellpack(A, beta, threshold);
            break;
        case ellblock:
            bml_add_identity_ellblock(A, beta, threshold);
            break;
        case ellsort:
            bml_add_identity_ellsort(A, beta, threshold);
            break;
        default:
            LOG_ERROR("unknown matrix type\n");
            break;
    }
}

/** Matrix addition.
 *
 * \f$ A \leftarrow \alpha A + \beta \mathrm{Id} \f$
 *
 * \ingroup add_group_C
 *
 * \param A Matrix A
 * \param alpha Scalar factor multiplied by A
 * \param beta Scalar factor multiplied by I
 * \param threshold Threshold for matrix addition
 */
void
bml_scale_add_identity(
    bml_matrix_t * A,
    double alpha,
    double beta,
    double threshold)
{
    switch (bml_get_type(A))
    {
        case dense:
            bml_scale_add_identity_dense(A, alpha, beta);
            break;
        case ellpack:
            bml_scale_add_identity_ellpack(A, alpha, beta, threshold);
            break;
        case ellblock:
            bml_scale_add_identity_ellblock(A, alpha, beta, threshold);
            break;
        case ellsort:
            bml_scale_add_identity_ellsort(A, alpha, beta, threshold);
            break;
        default:
            LOG_ERROR("unknown matrix type\n");
            break;
    }
}
