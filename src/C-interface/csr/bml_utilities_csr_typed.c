#include "../../macros.h"
#include "../../typed.h"
#include "../bml_logger.h"
#include "../bml_parallel.h"
#include "../bml_types.h"
#include "../bml_utilities.h"
#include "bml_types_csr.h"
#include "bml_setters_csr.h"
#include "bml_utilities_csr.h"

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
    bml_read_bml_matrix_csr) (
    bml_matrix_csr_t * A,
    char *filename)
{
    FILE *hFile;
    char header1[20], header2[20], header3[20], header4[20], header5[20];
    int hdimx, nnz, irow, icol;
    double real_part, imaginary_part;
    REAL_T value;

    hFile = fopen(filename, "r");

    // Read header
    if (fscanf(hFile, "%s %s %s %s %s", header1, header2, header3, header4,
               header5) != 5)
    {
        LOG_ERROR("read error on header\n");
    }

    LOG_DEBUG("Read: %s %s %s %s %s\n", header1, header2, header3, header4,
              header5);

    int symflag = strcmp(header5, "symmetric");

    // Read N, N, # of non-zeroes
    if (fscanf(hFile, "%d %d %d", &hdimx, &hdimx, &nnz) != 3)
    {
        LOG_ERROR("read error\n");
    }

    LOG_DEBUG("hdimx = %d, nnz = %d\n", hdimx, nnz);

    // Read in values
    for (int i = 0; i < nnz; i++)
    {
#if defined(SINGLE_REAL)
        if (fscanf(hFile, "%d %d %f\n", &irow, &icol, &value) != 3)
        {
            LOG_ERROR("read error\n");
        }
#elif defined(DOUBLE_REAL)
        if (fscanf(hFile, "%d %d %lf\n", &irow, &icol, &value) != 3)
        {
            LOG_ERROR("read error\n");
        }
#elif defined(SINGLE_COMPLEX)
        if (fscanf
            (hFile, "%d %d %lf %lf\n", &irow, &icol, &real_part,
             &imaginary_part) != 4)
        {
            LOG_ERROR("read error\n");
        }
        value = real_part + I * imaginary_part;
        LOG_DEBUG("read: %d %d %e %e %e\n", irow, icol, real_part,
                  imaginary_part, value);
#elif defined(DOUBLE_COMPLEX)
        if (fscanf
            (hFile, "%d %d %lf %lf\n", &irow, &icol, &real_part,
             &imaginary_part) != 4)
        {
            LOG_ERROR("read error\n");
        }
        value = real_part + I * imaginary_part;
#else
        LOG_ERROR("unknown precision\n");
#endif

        irow--;
        icol--;
        // set new element
        TYPED_FUNC(csr_set_row_element_new) (A->data_[irow], icol, &value);
        // Set symmetric value if necessary
        if (symflag == 0 && icol != irow)
        {
            TYPED_FUNC(csr_set_row_element_new) (A->data_[icol], irow,
                                                 &value);
        }
    }

    fclose(hFile);
}

/** Write a Matrix Market format file from a bml matrix.
 *
 *  \ingroup utilities_group
 *
 *  \param A The matrix to be written
 *  \param filename The Matrix Market format file
 */
void TYPED_FUNC(
    bml_write_bml_matrix_csr) (
    bml_matrix_csr_t * A,
    char *filename)
{
    FILE *mFile;
    int msum;

    int N = A->N_;

    REAL_T *vals = NULL;
    int *cols = NULL;

    // Only write from rank 0
    if (bml_printRank() != 1)
        return;

    mFile = fopen(filename, "w");

    // Write header
#if defined(SINGLE_REAL) || defined(DOUBLE_REAL)
    fprintf(mFile, "%%%%%%MatrixMarket matrix coordinate real general\n");
#elif defined(SINGLE_COMPLEX) || defined(DOUBLE_COMPLEX)
    fprintf(mFile, "%%%%%%MatrixMarket matrix coordinate complex general\n");
#endif

    // Collect number of non-zero elements
    // Write out matrix size as dense and number of non-zero elements
    msum = 0;
    for (int i = 0; i < N; i++)
    {
        msum += A->data_[i]->NNZ_;
    }
    fprintf(mFile, "%d %d %d\n", N, N, msum);

    // Write out non-zero elements
    for (int i = 0; i < N; i++)
    {
        cols = A->data_[i]->cols_;
        vals = (REAL_T *) A->data_[i]->vals_;
        const int annz = A->data_[i]->NNZ_;
        for (int pos = 0; pos < annz; pos++)
        {
#if defined(SINGLE_REAL)
            fprintf(mFile, "%d %d %20.15g\n", i + 1,
                    cols[pos] + 1, vals[pos]);
#elif defined(DOUBLE_REAL)
            fprintf(mFile, "%d %d %20.15lg\n", i + 1,
                    cols[pos] + 1, REAL_PART(vals[pos]));
#elif defined(SINGLE_COMPLEX)
            fprintf(mFile, "%d %d %20.15g %20.15g\n", i + 1,
                    cols[pos] + 1, REAL_PART(vals[pos]),
                    IMAGINARY_PART(vals[pos]));
#elif defined(DOUBLE_COMPLEX)
            fprintf(mFile, "%d %d %20.15lg %20.15lg\n", i + 1,
                    cols[pos] + 1, REAL_PART(vals[pos]),
                    IMAGINARY_PART(vals[pos]));
#else
            LOG_ERROR("unknown precision\n");
#endif
        }
    }

    fclose(mFile);
}
