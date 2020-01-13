#include "../../macros.h"
#include "../../typed.h"
#include "../bml_logger.h"
#include "../bml_parallel.h"
#include "../bml_types.h"
#include "../bml_utilities.h"
#include "bml_types_ellblock.h"
#include "bml_utilities_ellblock.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/** Read in a bml matrix from Matrix Market format.
 *
 *  \ingroup utilities_group
 *
 *  \param A The matrix to be read
 *  \param filename The Matrix Market format file
 */
void TYPED_FUNC(
    bml_read_bml_matrix_ellblock) (
    bml_matrix_ellblock_t * A,
    char *filename)
{
    assert(A->bsize[0] < 1e6);

    FILE *hFile;
    char header1[20], header2[20], header3[20], header4[20], header5[20];
    int hdimx, nnz, irow, icol;
    REAL_T val;

    int NB = A->NB;
    int MB = A->MB;
    REAL_T **A_ptr_value = (REAL_T **) A->ptr_value;
    int *A_indexb = A->indexb;
    int *A_nnzb = A->nnzb;
    int *A_bsize = A->bsize;

    hFile = fopen(filename, "r");

    // Read header
    if (fscanf(hFile, "%s %s %s %s %s", header1, header2, header3, header4,
               header5) != 5)
    {
        LOG_ERROR("read error on header\n");
    }

    int symflag = strcmp(header5, "symmetric");

    // Read N, N, # of non-zeroes
    if (fscanf(hFile, "%d %d %d", &hdimx, &hdimx, &nnz) != 3)
    {
        LOG_ERROR("read error\n");
    }

    char *FMT = "";
    switch (A->matrix_precision)
    {
        case single_real:
            FMT = "%d %d %g\n";
            break;
        case double_real:
            FMT = "%d %d %lg\n";
            break;
        case single_complex:
            FMT = "%d %d %g\n";
            break;
        case double_complex:
            FMT = "%d %d %lg\n";
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }

    // Read in values
    for (int i = 0; i < nnz; i++)
    {
        if (fscanf(hFile, FMT, &irow, &icol, &val) != 3)
        {
            LOG_ERROR("read error\n");
        }
        irow--;
        icol--;

        /*compute block indexes */
        int ib = 0;
        while (irow >= A_bsize[ib])
        {
            irow -= A_bsize[ib];
            ib++;
        }

        int jb = 0;
        while (icol >= A_bsize[jb])
        {
            icol -= A_bsize[jb];
            jb++;
        }

        int found_block = 0;
        for (int jp = 0; jp < A_nnzb[ib]; jp++)
        {
            int ind = ROWMAJOR(ib, jp, NB, MB);
            int jbb = A_indexb[ind];
            if (jbb == jb)
            {
                found_block = 1;
                REAL_T *A_value = A_ptr_value[ind];
                A_value[ROWMAJOR(irow, icol, A_bsize[ib], A_bsize[jb])] = val;
            }
        }
        //add block if needed
        if (found_block == 0)
        {
            //printf("Add new block %d, %d\n",ib, jb);
            assert(A_nnzb[ib] < MB);
            int nelements = A_bsize[ib] * A_bsize[jb];
            int ind = ROWMAJOR(ib, A_nnzb[ib], NB, MB);
            A_indexb[ind] = jb;
            A_nnzb[ib]++;
            A_ptr_value[ind] = calloc(nelements, sizeof(REAL_T));
            REAL_T *A_value = A_ptr_value[ind];
            assert(A_value != NULL);
            A_value[ROWMAJOR(irow, icol, A_bsize[ib], A_bsize[jb])] = val;
        }

        // Set symmetric value if necessary
        if (symflag == 0 && icol != irow)
        {
            int found_block = 0;
            for (int ip = 0; ip < A_nnzb[jb]; ip++)
            {
                int ind = ROWMAJOR(jb, ip, NB, MB);
                int ibb = A_indexb[ROWMAJOR(jb, ip, NB, MB)];
                if (ibb == ib)
                {
                    found_block = 1;
                    REAL_T *A_value = A_ptr_value[ind];
                    A_value[ROWMAJOR(icol, irow, A_bsize[jb], A_bsize[ib])] =
                        val;
                }
            }
            //add block if needed
            if (found_block == 0)
            {
                A_indexb[A_nnzb[jb]] = ib;
                int ind = ROWMAJOR(jb, A_nnzb[jb], NB, MB);
                A_nnzb[jb]++;
                REAL_T *A_value = A_ptr_value[ind];
                A_value[ROWMAJOR(icol, irow, A_bsize[jb], A_bsize[ib])] = val;
            }
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
    bml_write_bml_matrix_ellblock) (
    bml_matrix_ellblock_t * A,
    char *filename)
{
    FILE *mFile;
    int msum;

    int NB = A->NB;
    int MB = A->MB;

    int *A_indexb = A->indexb;
    int *A_nnzb = A->nnzb;
    int *bsize = A->bsize;
    REAL_T **A_ptr_value = (REAL_T **) A->ptr_value;

    // Only write from rank 0
    if (bml_printRank() != 1)
        return;

    mFile = fopen(filename, "w");

    // Write header
    fprintf(mFile, "%%%%%%MatrixMarket matrix coordinate real general\n");

    // Collect number of non-zero elements
    // Write out matrix size as dense and number of non-zero elements
    msum = 0;
    for (int ib = 0; ib < NB; ib++)
    {
        for (int jp = 0; jp < A_nnzb[ib]; jp++)
        {
            int kb = ROWMAJOR(ib, jp, NB, MB);
            int jb = A_indexb[kb];
            msum += bsize[ib] * bsize[jb];
        }
    }
    fprintf(mFile, "%d %d %d\n", A->N, A->N, msum);

    int *offset = malloc(NB * sizeof(int));
    offset[0] = 0;
    for (int ib = 1; ib < NB; ib++)
    {
        offset[ib] = offset[ib - 1] + bsize[ib - 1];
    }

    // Write out non-zero elements
    for (int ib = 0; ib < NB; ib++)
    {
        for (int jp = 0; jp < A_nnzb[ib]; jp++)
        {
            int kb = ROWMAJOR(ib, jp, NB, MB);
            int jb = A_indexb[kb];
            REAL_T *A_value = A_ptr_value[kb];
            for (int ii = 0; ii < bsize[ib]; ii++)
                for (int jj = 0; jj < bsize[jb]; jj++)
                {
                    fprintf(mFile, "%d %d %20.15e\n",
                            offset[ib] + 1 + ii,
                            offset[jb] + 1 + jj,
                            REAL_PART(A_value
                                      [ROWMAJOR
                                       (ii, jj, bsize[ib], bsize[jb])]));
                }
        }
    }
    free(offset);

    fclose(mFile);
}

double TYPED_FUNC(
    bml_norm_inf) (
    void *_v,
    int nrows,
    int ncols,
    int ld)
{
    REAL_T *v = _v;
    double norm = 0.;
    for (int i = 0; i < nrows; i++)
        for (int j = 0; j < ncols; j++)
        {
            double alpha = (double) ABS(v[ROWMAJOR(i, j, nrows, ld)]);
            if (alpha > norm)
                norm = alpha;
        }
    return norm;
}

double TYPED_FUNC(
    bml_norm_inf_fast) (
    void *_v,
    int n)
{
    REAL_T *v = _v;
    double norm = 0.;
    for (int i = 0; i < n; i++)
    {
        double alpha = (double) ABS(v[i]);
        if (alpha > norm)
            norm = alpha;
    }
    return norm;
}

double TYPED_FUNC(
    bml_sum_squares) (
    void *_v,
    int nrows,
    int ncols,
    int ld)
{
    REAL_T *v = _v;
    double n2 = 0.;
    for (int i = 0; i < nrows; i++)
        for (int j = 0; j < ncols; j++)
        {
            double alpha = (double) ABS(v[ROWMAJOR(i, j, nrows, ld)]);
            n2 = n2 + alpha * alpha;
        }
    return n2;
}
