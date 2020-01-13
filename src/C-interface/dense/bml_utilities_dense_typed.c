#ifdef BML_USE_MAGMA
#include "magma_v2.h"
#endif

#include "../../typed.h"
#include "../../macros.h"
#include "../bml_logger.h"
#include "../bml_utilities.h"
#include "../bml_allocate.h"
#include "bml_types_dense.h"
#include "bml_utilities_dense.h"
#include "bml_export_dense.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/** Read in a bml matrix from Matrix Market format.
 *
 *  \ingroup utilities_group
 *
 *  \param A The matrix to be read
 *  \param filename The Matrix Market format file
 */
void TYPED_FUNC(
    bml_read_bml_matrix_dense) (
    bml_matrix_dense_t * A,
    char *filename)
{
    FILE *matrix_file;
    char header1[20], header2[20], header3[20], header4[20], header5[20];
    int hdimx, nnz, irow, icol;
    int values_read;
    double real_part, imaginary_part;

    int N = A->N;

    REAL_T *A_value = A->matrix;

    matrix_file = fopen(filename, "r");

    // Read header
    if ((values_read = fscanf
         (matrix_file, "%s %s %s %s %s", header1, header2, header3, header4,
          header5)) < 5)
    {
        LOG_ERROR("Line 1, expected 5 entries, read %d\n", values_read);
    }

    // Read N, N, # of non-zeroes
    if ((values_read =
         fscanf(matrix_file, "%d %d %d", &hdimx, &hdimx, &nnz)) < 3)
    {
        LOG_ERROR("Line 2, Expected 3 entries, read %d\n", values_read);
    }

    // Read in values
    for (int i = 0; i < nnz; i++)
    {
        switch (A->matrix_precision)
        {
            case single_real:
                if ((values_read = fscanf
                     (matrix_file, "%d %d %lg\n", &irow, &icol,
                      &real_part)) < 3)
                {
                    LOG_ERROR("Line %d, expected 3 entries, read %d\n", i + 3,
                              values_read);
                }
                A_value[ROWMAJOR(irow - 1, icol - 1, N, N)] = real_part;
                break;
            case double_real:
                if ((values_read = fscanf
                     (matrix_file, "%d %d %lg\n", &irow, &icol,
                      &real_part)) < 3)
                {
                    LOG_ERROR("Line %d, expected 3 entries, read %d\n", i + 3,
                              values_read);
                }
                A_value[ROWMAJOR(irow - 1, icol - 1, N, N)] = real_part;
                break;
            case single_complex:
                if ((values_read = fscanf
                     (matrix_file, "%d %d %lg %lg\n", &irow, &icol,
                      &real_part, &imaginary_part)) < 4)
                {
                    LOG_ERROR("Line %d, expected 4 entries, read %d\n", i + 3,
                              values_read);
                }
                A_value[ROWMAJOR(irow - 1, icol - 1, N, N)] =
                    real_part + I * imaginary_part;
                break;
            case double_complex:
                if ((values_read = fscanf
                     (matrix_file, "%d %d %lg %lg\n", &irow, &icol,
                      &real_part, &imaginary_part)) < 4)
                {
                    LOG_ERROR("Line %d, expected 4 entries, read %d\n", i + 3,
                              values_read);
                }
                A_value[ROWMAJOR(irow - 1, icol - 1, N, N)] =
                    real_part + I * imaginary_part;
                break;
            default:
                LOG_ERROR("unknown precision\n");
                break;
        }
    }

    fclose(matrix_file);
}

/** Write a Matrix Market format file from a bml matrix.
 *
 *  \ingroup utilities_group
 *
 * Note, all matrix elements are written, even the ones that are zero.
 *
 * \param A The matrix to be written
 * \param filename The Matrix Market format file
 */
void TYPED_FUNC(
    bml_write_bml_matrix_dense) (
    bml_matrix_dense_t * A,
    char *filename)
{
    FILE *matrix_file;

    int N = A->N;
    int msum = N * N;

    REAL_T *A_value = A->matrix;

    matrix_file = fopen(filename, "w");

    // Write header
    if (A->matrix_precision == single_complex
        || A->matrix_precision == double_complex)
    {
        fprintf(matrix_file,
                "%%%%%%MatrixMarket matrix coordinate complex general\n");
    }
    else
    {
        fprintf(matrix_file,
                "%%%%%%MatrixMarket matrix coordinate real general\n");
    }

    // Write out matrix size as dense and number of non-zero elements
    fprintf(matrix_file, "%d %d %d\n", N, N, msum);

    // Write out non-zero elements
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            switch (A->matrix_precision)
            {
                case single_real:
                    fprintf(matrix_file, "%d %d %20.15g\n", i + 1, j + 1,
                            REAL_PART(A_value[ROWMAJOR(i, j, N, N)]));
                    break;
                case double_real:
                    fprintf(matrix_file, "%d %d %20.15lg\n", i + 1, j + 1,
                            REAL_PART(A_value[ROWMAJOR(i, j, N, N)]));
                    break;
                case single_complex:
                    fprintf(matrix_file, "%d %d %20.15g %20.15g\n", i + 1,
                            j + 1, REAL_PART(A_value[ROWMAJOR(i, j, N, N)]),
                            IMAGINARY_PART(A_value[ROWMAJOR(i, j, N, N)]));
                    break;
                case double_complex:
                    fprintf(matrix_file, "%d %d %20.15lg %20.15lg\n", i + 1,
                            j + 1, REAL_PART(A_value[ROWMAJOR(i, j, N, N)]),
                            IMAGINARY_PART(A_value[ROWMAJOR(i, j, N, N)]));
                    break;
                default:
                    LOG_ERROR("unknown precision\n");
                    break;
            }
        }
    }
    fclose(matrix_file);
}

void TYPED_FUNC(
    bml_print_bml_matrix_dense) (
    bml_matrix_dense_t * A,
    int i_l,
    int i_u,
    int j_l,
    int j_u)
{
#ifdef BML_USE_MAGMA
    //copy matrix data from GPU to CPU so that we can use
    //bml_print_dense_matrix function with lower/upper row/column
    //index
    REAL_T *A_matrix = bml_allocate_memory(sizeof(REAL_T) * A->N * A->N);
    MAGMA(getmatrix) (A->N, A->N,
                      A->matrix, A->ld, (MAGMA_T *) A_matrix, A->N, A->queue);
#else
    REAL_T *A_matrix = (REAL_T *) A->matrix;
#endif

    bml_print_dense_matrix(A->N, A->matrix_precision, dense_row_major,
                           A_matrix, i_l, i_u, j_l, j_u);

#ifdef BML_USE_MAGMA
    bml_free_memory(A_matrix);
#endif
}
