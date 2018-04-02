#include "../../macros.h"
#include "../dense/bml_types_dense.h"
#include "bml_logger.h"
#include "bml_submatrix.h"
#include "bml_submatrix_ellpack.h"
#include "bml_types.h"
#include "bml_types_ellpack.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

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
 * \param vsize Size of core_halo_index and core_pos
 * \param double_jump_flag Flag to use double jump (0=no, 1=yes)
 */
void
bml_matrix2submatrix_index_ellpack(
    const bml_matrix_ellpack_t * A,
    const bml_matrix_ellpack_t * B,
    const int *nodelist,
    const int nsize,
    int *core_halo_index,
    int *vsize,
    const int double_jump_flag)
{
    switch (A->matrix_precision)
    {
        case single_real:
            bml_matrix2submatrix_index_ellpack_single_real(A, B, nodelist,
                                                           nsize,
                                                           core_halo_index,
                                                           vsize,
                                                           double_jump_flag);
            break;
        case double_real:
            bml_matrix2submatrix_index_ellpack_double_real(A, B, nodelist,
                                                           nsize,
                                                           core_halo_index,
                                                           vsize,
                                                           double_jump_flag);
            break;
        case single_complex:
            bml_matrix2submatrix_index_ellpack_single_complex(A, B, nodelist,
                                                              nsize,
                                                              core_halo_index,
                                                              vsize,
                                                              double_jump_flag);
            break;
        case double_complex:
            bml_matrix2submatrix_index_ellpack_double_complex(A, B, nodelist,
                                                              nsize,
                                                              core_halo_index,
                                                              vsize,
                                                              double_jump_flag);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
}

/** Determine element indices for submatrix, given a set of nodes/orbitals.
 *
 * \ingroup submatrix_group_C
 *
 * \param B Graph matrix B
 * \param nodelist List of node/orbital indeces
 * \param nsize Size of nodelist
 * \param core_halo_index List of core+halo indeces
 * \param vsize Size of core_halo_index and core_pos
 * \param double_jump_flag Flag to use double jump (0=no, 1=yes)
 */
void
bml_matrix2submatrix_index_graph_ellpack(
    const bml_matrix_ellpack_t * B,
    const int *nodelist,
    const int nsize,
    int *core_halo_index,
    int *vsize,
    const int double_jump_flag)
{
    switch (B->matrix_precision)
    {
        case single_real:
            bml_matrix2submatrix_index_graph_ellpack_single_real(B, nodelist,
                                                                 nsize,
                                                                 core_halo_index,
                                                                 vsize,
                                                                 double_jump_flag);
            break;
        case double_real:
            bml_matrix2submatrix_index_graph_ellpack_single_real(B, nodelist,
                                                                 nsize,
                                                                 core_halo_index,
                                                                 vsize,
                                                                 double_jump_flag);
            break;
        case single_complex:
            bml_matrix2submatrix_index_graph_ellpack_double_complex(B,
                                                                    nodelist,
                                                                    nsize,
                                                                    core_halo_index,
                                                                    vsize,
                                                                    double_jump_flag);
            break;
        case double_complex:
            bml_matrix2submatrix_index_graph_ellpack_double_complex(B,
                                                                    nodelist,
                                                                    nsize,
                                                                    core_halo_index,
                                                                    vsize,
                                                                    double_jump_flag);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
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
void
bml_matrix2submatrix_ellpack(
    const bml_matrix_ellpack_t * A,
    bml_matrix_dense_t * B,
    const int *core_halo_index,
    const int lsize)
{
    switch (A->matrix_precision)
    {
        case single_real:
            bml_matrix2submatrix_ellpack_single_real(A, B, core_halo_index,
                                                     lsize);
            break;
        case double_real:
            bml_matrix2submatrix_ellpack_double_real(A, B, core_halo_index,
                                                     lsize);
            break;
        case single_complex:
            bml_matrix2submatrix_ellpack_single_complex(A, B, core_halo_index,
                                                        lsize);
            break;
        case double_complex:
            bml_matrix2submatrix_ellpack_double_complex(A, B, core_halo_index,
                                                        lsize);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
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
 * \param threshold Threshold for elements
 */
void
bml_submatrix2matrix_ellpack(
    const bml_matrix_dense_t * A,
    bml_matrix_ellpack_t * B,
    const int *core_halo_index,
    const int lsize,
    const int llsize,
    const double threshold)
{
    switch (A->matrix_precision)
    {
        case single_real:
            bml_submatrix2matrix_ellpack_single_real(A, B, core_halo_index,
                                                     lsize, llsize,
                                                     threshold);
            break;
        case double_real:
            bml_submatrix2matrix_ellpack_double_real(A, B, core_halo_index,
                                                     lsize, llsize,
                                                     threshold);
            break;
        case single_complex:
            bml_submatrix2matrix_ellpack_single_complex(A, B, core_halo_index,
                                                        lsize,
                                                        llsize, threshold);
            break;
        case double_complex:
            bml_submatrix2matrix_ellpack_double_complex(A, B, core_halo_index,
                                                        lsize,
                                                        llsize, threshold);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
}

/** Get vector from matrix.
 *
 * \ingroup submatrix_group_C
 *
 * \param A Matrix A
 * \param jj Index set
 * \param irow Which row
 * \param colCnt Number of columns
 * \param rvalue Returned vector
 */
void *
bml_getVector_ellpack(
    const bml_matrix_ellpack_t * A,
    const int *jj,
    const int irow,
    const int colCnt)
{
    switch (A->matrix_precision)
    {
        case single_real:
            return bml_getVector_ellpack_single_real(A, jj, irow, colCnt);
            break;
        case double_real:
            return bml_getVector_ellpack_double_real(A, jj, irow, colCnt);
            break;
        case single_complex:
            return bml_getVector_ellpack_single_complex(A, jj, irow, colCnt);
            break;
        case double_complex:
            return bml_getVector_ellpack_double_complex(A, jj, irow, colCnt);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
    return NULL;
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
bml_matrix_ellpack_t *
bml_group_matrix_ellpack(
    const bml_matrix_ellpack_t * A,
    const int *hindex,
    const int ngroups,
    const double threshold)
{
    switch (A->matrix_precision)
    {
        case single_real:
            return bml_group_matrix_ellpack_single_real(A, hindex, ngroups,
                                                        threshold);
            break;
        case double_real:
            return bml_group_matrix_ellpack_double_real(A, hindex, ngroups,
                                                        threshold);
            break;
        case single_complex:
            return bml_group_matrix_ellpack_single_complex(A, hindex, ngroups,
                                                           threshold);
            break;
        case double_complex:
            return bml_group_matrix_ellpack_double_complex(A, hindex, ngroups,
                                                           threshold);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
    return NULL;
}

int
sortById(
    const void *a,
    const void *b)
{
    int aId = *((int *) a);
    int bId = *((int *) b);

    if (aId < bId)
        return -1;
    else if (aId == bId)
        return 0;
    else
        return 1;
}

/** Assemble adjacency structure from matrix.
 *
 * \ingroup submatrix_group_C
 *
 * \param A Matrix A
 * \param xadj Index of each row in adjncy
 * \param adjncy Adjacency vector
 * \param base_flag Return 0- or 1-based
 */
void
bml_adjacency_ellpack(
    const bml_matrix_ellpack_t * A,
    int *xadj,
    int *adjncy,
    const int base_flag)
{
    int A_N = A->N;
    int A_M = A->M;

    int *A_nnz = A->nnz;
    int *A_index = A->index;

    int j;
    int check;

    xadj[0] = 0;

    // Check if diagonal elements are included
    check = 0;
    for (int i = 0; i < A_nnz[0]; i++)
    {
        if (A_index[ROWMAJOR(0, i, A_N, A_M)] == 0)
        {
            check = 1;
            break;
        }
    }

    for (int i = 1; i < A_N + 1; i++)
    {
        if (check == 1)
            xadj[i] = xadj[i - 1] + A_nnz[i - 1] - 1;
        else
            xadj[i] = xadj[i - 1] + A_nnz[i - 1];
    }

#pragma omp parallel for default(none) \
    private(j) \
    shared(A_N, A_M, A_index, A_nnz, xadj, adjncy)
    for (int i = 0; i < A_N; i++)
    {
        j = xadj[i];
        for (int jj = 0; jj < A_nnz[i]; jj++)
        {
            if (A_index[ROWMAJOR(i, jj, A_N, A_M)] != i)
            {
                adjncy[j] = A_index[ROWMAJOR(i, jj, A_N, A_M)];
                j++;
            }
        }
        //assert(j == xadj[i+1]);
    }

#pragma omp parallel for default(none) \
    shared(A_N, xadj, adjncy)
    for (int i = 0; i < A_N; i++)
    {
        qsort(&adjncy[xadj[i]], xadj[i + 1] - xadj[i], sizeof(int), sortById);
    }

    // Add 1 for 1-based
    if (base_flag == 1)
    {
#pragma omp parallel for default(none) \
    shared(xadj, A_N, adjncy)
        for (int i = 0; i < A_N; i++)
        {
            for (int j = xadj[i]; j < xadj[i + 1]; j++)
            {
                adjncy[j] += 1;
            }
        }
#pragma omp parallel for default(none) \
    shared(xadj, A_N)
        for (int i = 0; i < A_N + 1; i++)
        {
            xadj[i] += 1;
        }
    }
}

/** Assemble adjacency structure from matrix based on groups of rows.
 *
 * \ingroup submatrix_group_C
 *
 * \param A Matrix A
 * \param hindex Indeces of nodes
 * \param nnodes Number of groups
 * \param xadj Index of each row in adjncy
 * \param adjncy Adjacency vector
 * \param base_flag Return 0- or 1-based
 */
void
bml_adjacency_group_ellpack(
    const bml_matrix_ellpack_t * A,
    const int *hindex,
    const int nnodes,
    int *xadj,
    int *adjncy,
    const int base_flag)
{
    int A_N = A->N;
    int A_M = A->M;

    int *A_nnz = A->nnz;
    int *A_index = A->index;

    int *hnode = malloc(nnodes * sizeof(int));
    for (int i = 0; i < nnodes; i++)
    {
        hnode[i] = hindex[i] - 1;
    }

    // Determine number of adjacent atoms per atom
    xadj[0] = 0;
    for (int i = 1; i < nnodes + 1; i++)
    {
        int hcount = 0;
        for (int j = 0; j < nnodes; j++)
        {
            for (int k = 0; k < A_nnz[hnode[i - 1]]; k++)
            {
                if (hnode[j] == A_index[ROWMAJOR(hnode[i - 1], k, A_N, A_M)])
                {
                    hcount++;
                    break;
                }
            }
        }

        xadj[i] = xadj[i - 1] + hcount;
    }

    // Fill in adjacent atoms
#pragma omp parallel for default(none) \
    shared(A_N, A_M, A_index, A_nnz) \
    shared(xadj, adjncy, hnode)
    for (int i = 0; i < nnodes; i++)
    {
        int ll = xadj[i];

        for (int j = 0; j < nnodes; j++)
        {
            for (int k = 0; k < A_nnz[hnode[i]]; k++)
            {
                if (hnode[j] == A_index[ROWMAJOR(hnode[i], k, A_N, A_M)])
                {
                    //adjncy[ll] = hnode[j];
                    adjncy[ll] = j;
                    ll++;
                    break;
                }
            }
        }
    }

    // Add 1 for 1-based
    if (base_flag == 1)
    {
#pragma omp parallel for default(none) \
    shared(xadj, A_N, adjncy)
        for (int i = 0; i <= xadj[nnodes]; i++)
        {
            adjncy[i] += 1;
        }
#pragma omp parallel for default(none) \
    shared(xadj, A_N)
        for (int i = 0; i < nnodes + 1; i++)
        {
            xadj[i] += 1;
        }
    }

}
