#include "../../macros.h"
#include "../../typed.h"
#include "../dense/bml_allocate_dense.h"
#include "bml_allocate.h"
#include "bml_allocate_ellsort.h"
#include "bml_logger.h"
#include "bml_submatrix.h"
#include "bml_submatrix_ellsort.h"
#include "bml_types.h"
#include "bml_types_ellsort.h"

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/** Determine element indices for submatrix, given a set of nodes/orbitals.
 *
 * \ingroup submatrix_group_C
 *
 * \param A Hamiltonian matrix A
 * \param B Graph matrix B
 * \param nodelist List of node/orbital indeces
 * \param nsize Size of nodelist
 * \param core_halo_index List of core+halo indeces
 * \param vsize Size of core_halo_index and number of cores
 * \param double_jump_flag Flag to use double jump (0=no, 1=yes)
 */
void TYPED_FUNC(
    bml_matrix2submatrix_index_ellsort) (
    const bml_matrix_ellsort_t * A,
    const bml_matrix_ellsort_t * B,
    const int *nodelist,
    const int nsize,
    int *core_halo_index,
    int *vsize,
    const int double_jump_flag)
{
    int l, ll, ii, ls, k;
    int A_N = A->N;
    int A_M = A->M;
    int *A_nnz = A->nnz;
    int *A_index = A->index;
    int B_N = B->N;
    int B_M = B->M;
    int *B_nnz = B->nnz;
    int *B_index = B->index;

    int ix[A_N];

    memset(ix, 0, A_N * sizeof(int));

    l = 0;
    ll = 0;

    // Cores are first followed by halos
    for (int j = 0; j < nsize; j++)
    {
        ii = nodelist[j];
        if (ix[ii] == 0)
        {
            ix[ii] = ii + 1;
            core_halo_index[l] = ii;
            l++;
            ll++;
        }

    }

    // Collect halo indeces from graph
    for (int j = 0; j < nsize; j++)
    {
        ii = nodelist[j];

        for (int jp = 0; jp < B_nnz[ii]; jp++)
        {
            k = B_index[ROWMAJOR(ii, jp, B_N, B_M)];
            if (ix[k] == 0)
            {
                ix[k] = ii + 1;
                core_halo_index[l] = k;
                l++;
            }
        }
    }

    // Add more halo elements from H
    for (int j = 0; j < nsize; j++)
    {
        ii = nodelist[j];

        for (int jp = 0; jp < A_nnz[ii]; jp++)
        {
            k = A_index[ROWMAJOR(ii, jp, A_N, A_M)];
            if (ix[k] == 0)
            {
                ix[k] = ii + 1;
                core_halo_index[l] = k;
                l++;
            }
        }
    }

    // Perform a "double jump" for extra halo elements
    // based on graph, like performing a symbolic X^2
    if (double_jump_flag == 1)
    {
        ls = l;
        for (int j = 0; j < ls; j++)
        {
            ii = core_halo_index[j];

            for (int jp = 0; jp < B_nnz[ii]; jp++)
            {
                k = B_index[ROWMAJOR(ii, jp, B_N, B_M)];
                if (ix[k] == 0)
                {
                    ix[k] = ii + 1;
                    core_halo_index[l] = k;
                    l++;
                }
            }
        }
    }

    vsize[0] = l;
    vsize[1] = ll;
}

/** Determine element indices for submatrix, given a set of nodes/orbitals.
 *
 * \ingroup submatrix_group_C
 *
 * \param B Graph matrix B
 * \param nodelist List of node/orbital indeces
 * \param nsize Size of nodelist
 * \param core_halo_index List of core+halo indeces
 * \param vsize Size of core_halo_index and number of cores
 * \param double_jump_flag Flag to use double jump (0=no, 1=yes)
 */
void TYPED_FUNC(
    bml_matrix2submatrix_index_graph_ellsort) (
    const bml_matrix_ellsort_t * B,
    const int *nodelist,
    const int nsize,
    int *core_halo_index,
    int *vsize,
    const int double_jump_flag)
{
    int l, ll, ii, ls, k;
    int B_N = B->N;
    int B_M = B->M;
    int *B_index = B->index;
    int *B_nnz = B->nnz;
    int ix[B_N];

    memset(ix, 0, B_N * sizeof(int));

    l = 0;
    ll = 0;

    // Cores are first followed by halos
    for (int j = 0; j < nsize; j++)
    {
        ii = nodelist[j];
        if (ix[ii] == 0)
        {
            ix[ii] = ii + 1;
            core_halo_index[l] = ii;
            l++;
            ll++;
        }
    }

    // Collext halo indeces from graph
    for (int j = 0; j < nsize; j++)
    {
        ii = nodelist[j];

        for (int jp = 0; jp < B_nnz[ii]; jp++)
        {
            k = B_index[ROWMAJOR(ii, jp, B_N, B_M)];
            if (ix[k] == 0)
            {
                ix[k] = ii + 1;
                core_halo_index[l] = k;
                l++;
            }
        }
    }

    // Use graph for double jumps
    if (double_jump_flag == 1)
    {
        ls = l;
        for (int j = 0; j < ls; j++)
        {
            ii = core_halo_index[j];

            for (int jp = 0; jp < B_nnz[ii]; jp++)
            {
                k = B_index[ROWMAJOR(ii, jp, B_N, B_M)];
                if (ix[k] == 0)
                {
                    ix[k] = ii + 1;
                    core_halo_index[l] = k;
                    l++;
                }
            }
        }
    }

    vsize[0] = l;
    vsize[1] = ll;
}

/** Extract a submatrix from a matrix given a set of core+halo rows.
 *
 * \ingroup submatrix_group_C
 *
 * \param A Matrix A
 * \param B Submatrix B
 * \param core_halo_index Set of row indeces for submatrix
 * \param llsize Number of indeces
 */
void TYPED_FUNC(
    bml_matrix2submatrix_ellsort) (
    const bml_matrix_ellsort_t * A,
    bml_matrix_dense_t * B,
    const int *core_halo_index,
    const int lsize)
{
    REAL_T *rvalue;

    int B_N = B->N;
    REAL_T *B_matrix = B->matrix;

#pragma omp parallel for \
    default(none) \
    private(rvalue) \
    shared(core_halo_index) \
    shared(A, B_matrix, B_N)
    for (int jb = 0; jb < lsize; jb++)
    {
        rvalue = TYPED_FUNC(bml_getVector_ellsort) (A, core_halo_index,
                                                    core_halo_index[jb],
                                                    lsize);
        for (int j = 0; j < lsize; j++)
        {
            B_matrix[ROWMAJOR(jb, j, B_N, B_N)] = rvalue[j];
        }

        free(rvalue);
    }
}

/** Assemble submatrix into a full matrix based on core+halo indeces.
 *
 * \ingroup submatrix_group_C
 *
 * \param A Submatrix A
 * \param B Matrix B
 * \param core_halo_index Set of submatrix row indeces
 * \param lsize Number of indeces
 * \param llsize Number of core positions
 */
void TYPED_FUNC(
    bml_submatrix2matrix_ellsort) (
    const bml_matrix_dense_t * A,
    bml_matrix_ellsort_t * B,
    const int *core_halo_index,
    const int lsize,
    const int llsize,
    const double threshold)
{
    int A_N = A->N;
    REAL_T *A_matrix = A->matrix;

    int B_N = B->N;
    int B_M = B->M;
    int *B_nnz = B->nnz;
    int *B_index = B->index;
    REAL_T *B_value = B->value;

    int ii, icol;

#pragma omp parallel for \
    default(none) \
    private(ii, icol) \
    shared(core_halo_index) \
    shared(A_N, A_matrix) \
    shared(B_N, B_M, B_nnz, B_index, B_value)
    for (int ja = 0; ja < llsize; ja++)
    {
        ii = core_halo_index[ja];

        icol = 0;
        for (int jb = 0; jb < lsize; jb++)
        {
            if (ABS(A_matrix[ROWMAJOR(ja, jb, A_N, A_N)]) > threshold)
            {
                B_index[ROWMAJOR(ii, icol, B_N, B_M)] = core_halo_index[jb];
                B_value[ROWMAJOR(ii, icol, B_N, B_M)] =
                    A_matrix[ROWMAJOR(ja, jb, A_N, A_N)];
                icol++;
            }
        }
        if (icol > B_M)
        {
            LOG_ERROR("Number of non-zeroes per row >= M, Increase M\n");
        }

        B_nnz[ii] = icol;
    }
}

// Get matching vector of values
void *TYPED_FUNC(
    bml_getVector_ellsort) (
    const bml_matrix_ellsort_t * A,
    const int *jj,
    const int irow,
    const int colCnt)
{
    const REAL_T ZERO = 0.0;

    int A_N = A->N;
    int A_M = A->M;
    int *A_nnz = A->nnz;
    int *A_index = A->index;
    REAL_T *A_value = A->value;
    REAL_T *rvalue = malloc(colCnt * sizeof(REAL_T));

    for (int i = 0; i < colCnt; i++)
    {
        for (int j = 0; j < A_nnz[irow]; j++)
        {
            if (A_index[ROWMAJOR(irow, j, A_N, A_M)] == jj[i])
            {
                rvalue[i] = A_value[ROWMAJOR(irow, j, A_N, A_M)];
                break;
            }
            rvalue[i] = ZERO;
        }
    }
    return rvalue;
}

/** Assemble matrix based on groups of rows from a matrix.
 *
 * \ingroup submatrix_group_C
 *
 * \param A Matrix A
 * \param hindex Indeces of nodes
 * \param ngroups Number of groups
 * \param threshold Threshold for graph
 */
bml_matrix_ellsort_t *TYPED_FUNC(
    bml_group_matrix_ellsort) (
    const bml_matrix_ellsort_t * A,
    const int *hindex,
    const int ngroups,
    const double threshold)
{
    int A_N = A->N;
    int A_M = A->M;
    int *A_index = A->index;
    int *A_nnz = A->nnz;
    REAL_T *A_value = A->value;

    int ix[ngroups];
    int hnode[A_N];
    int hend;

    bml_matrix_ellsort_t *B =
        TYPED_FUNC(bml_noinit_matrix_ellsort) (ngroups, ngroups,
                                               A->distribution_mode);

    int B_N = B->N;
    int B_M = B->M;
    int *B_index = B->index;
    int *B_nnz = B->nnz;
    REAL_T *B_value = B->value;

#pragma omp parallel for default(none) \
    private(hend) \
    shared(hindex, hnode, A_N)
    for (int i = 0; i < ngroups; i++)
    {
        hend = hindex[i + 1] - 1;
        if (i == ngroups - 1)
            hend = A_N;
        for (int j = hindex[i] - 1; j < hend; j++)
        {
            hnode[j] = i;
        }
    }

#pragma omp parallel for \
    default(none) \
    private(ix, hend) \
    shared(hindex, hnode) \
    shared(A_nnz, A_index, A_value, A_N, A_M) \
    shared(B_nnz, B_index, B_value, B_N, B_M)
    for (int i = 0; i < B_N; i++)
    {
        memset(ix, 0, sizeof(int) * ngroups);
        B_nnz[i] = 0;
        hend = hindex[i + 1] - 1;
        if (i == B_N - 1)
            hend = A_N;
        for (int j = hindex[i] - 1; j < hend; j++)
        {
            for (int k = 0; k < A_nnz[j]; k++)
            {
                int ii = hnode[A_index[ROWMAJOR(j, k, A_N, A_M)]];
                if (ix[ii] == 0 &&
                    is_above_threshold(A_value[ROWMAJOR(j, k, A_N, A_M)],
                                       threshold))
                {
                    ix[ii] = i + 1;
                    B_index[ROWMAJOR(i, B_nnz[i], B_N, B_M)] = ii;
                    B_value[ROWMAJOR(i, B_nnz[i], B_N, B_M)] = 1.0;
                    B_nnz[i]++;
                }
            }
        }
    }

    return B;
}
