/** \file */

#ifndef __MACROS_H
#define __MACROS_H

/** Row major access.
 *
 * \param i Row index.
 * \param j Column index.
 * \param M Number of rows.
 * \param N Number of columns.
 */
#define ROWMAJOR(i, j, M, N) (i) * (N) + (j)

/** Column major access.
 *
 * \param i Row index.
 * \param j Column index.
 * \param M Number of rows.
 * \param N Number of columns.
 */
#define COLMAJOR(i, j, M, N) (i) + (j) * (M)

/** The linear index into the dense matrix blocks of the BML.
 *
 * \param i Row index.
 * \param j Column index.
 * \param M Number of rows.
 * \param N Number of columns.
 */
#define LINEARINDEX(i, j, M, N) ROWMAJOR(i, j, M, N)

/** The min of two arguments.
 *
 * \param a The first argument.
 * \param b The second argument.
 * \return The smaller of the two input arguments.
 */
#define MIN(a, b) ((a) < (b) ? a : b)

/** The max of two arguments.
 *
 * \param a The first argument.
 * \param b The second argument.
 * \return The larger of the two input arguments.
 */
#define MAX(a, b) ((a) > (b) ? a : b)

#endif
