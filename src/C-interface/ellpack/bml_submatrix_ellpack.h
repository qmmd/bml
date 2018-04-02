#ifndef __BML_SUBMATRIX_ELLPACK_H
#define __BML_SUBMATRIX_ELLPACK_H

#include "bml_types_ellpack.h"
#include "../dense/bml_types_dense.h"

#include <complex.h>

void bml_matrix2submatrix_index_ellpack(
    const bml_matrix_ellpack_t * A,
    const bml_matrix_ellpack_t * B,
    const int *nodelist,
    const int nsize,
    int *core_halo_index,
    int *vsize,
    const int double_jump_flag);

void bml_matrix2submatrix_index_ellpack_single_real(
    const bml_matrix_ellpack_t * A,
    const bml_matrix_ellpack_t * B,
    const int *nodelist,
    const int nsize,
    int *core_halo_index,
    int *vsize,
    const int double_jump_flag);

void bml_matrix2submatrix_index_ellpack_double_real(
    const bml_matrix_ellpack_t * A,
    const bml_matrix_ellpack_t * B,
    const int *nodelist,
    const int nsize,
    int *core_halo_index,
    int *vsize,
    const int double_jump_flag);

void bml_matrix2submatrix_index_ellpack_single_complex(
    const bml_matrix_ellpack_t * A,
    const bml_matrix_ellpack_t * B,
    const int *nodelist,
    const int nsize,
    int *core_halo_index,
    int *vsize,
    const int double_jump_flag);

void bml_matrix2submatrix_index_ellpack_double_complex(
    const bml_matrix_ellpack_t * A,
    const bml_matrix_ellpack_t * B,
    const int *nodelist,
    const int nsize,
    int *core_halo_index,
    int *vsize,
    const int double_jump_flag);

void bml_matrix2submatrix_index_graph_ellpack(
    const bml_matrix_ellpack_t * B,
    const int *nodelist,
    const int nsize,
    int *core_halo_index,
    int *vsize,
    const int double_jump_flag);

void bml_matrix2submatrix_index_graph_ellpack_single_real(
    const bml_matrix_ellpack_t * B,
    const int *nodelist,
    const int nsize,
    int *core_halo_index,
    int *vsize,
    const int double_jump_flag);

void bml_matrix2submatrix_index_graph_ellpack_double_real(
    const bml_matrix_ellpack_t * B,
    const int *nodelist,
    const int nsize,
    int *core_halo_index,
    int *vsize,
    const int double_jump_flag);

void bml_matrix2submatrix_index_graph_ellpack_single_complex(
    const bml_matrix_ellpack_t * B,
    const int *nodelist,
    const int nsize,
    int *core_halo_index,
    int *vsize,
    const int double_jump_flag);

void bml_matrix2submatrix_index_graph_ellpack_double_complex(
    const bml_matrix_ellpack_t * B,
    const int *nodelist,
    const int nsize,
    int *core_halo_index,
    int *vsize,
    const int double_jump_flag);

void bml_matrix2submatrix_ellpack(
    const bml_matrix_ellpack_t * A,
    bml_matrix_dense_t * B,
    const int *core_halo_index,
    const int lsize);

void bml_matrix2submatrix_ellpack_single_real(
    const bml_matrix_ellpack_t * A,
    bml_matrix_dense_t * B,
    const int *core_halo_index,
    const int lsize);

void bml_matrix2submatrix_ellpack_double_real(
    const bml_matrix_ellpack_t * A,
    bml_matrix_dense_t * B,
    const int *core_halo_index,
    const int lsize);

void bml_matrix2submatrix_ellpack_single_complex(
    const bml_matrix_ellpack_t * A,
    bml_matrix_dense_t * B,
    const int *core_halo_index,
    const int lsize);

void bml_matrix2submatrix_ellpack_double_complex(
    const bml_matrix_ellpack_t * A,
    bml_matrix_dense_t * B,
    const int *core_halo_index,
    const int lsize);

void bml_submatrix2matrix_ellpack(
    const bml_matrix_dense_t * A,
    bml_matrix_ellpack_t * B,
    const int *core_halo_index,
    const int lsize,
    const int llsize,
    const double threshold);

void bml_submatrix2matrix_ellpack_single_real(
    const bml_matrix_dense_t * A,
    bml_matrix_ellpack_t * B,
    const int *core_halo_index,
    const int lsize,
    const int llsize,
    const double threshold);

void bml_submatrix2matrix_ellpack_double_real(
    const bml_matrix_dense_t * A,
    bml_matrix_ellpack_t * B,
    const int *core_halo_index,
    const int lsize,
    const int llsize,
    const double threshold);

void bml_submatrix2matrix_ellpack_single_complex(
    const bml_matrix_dense_t * A,
    bml_matrix_ellpack_t * B,
    const int *core_halo_index,
    const int lsize,
    const int llsize,
    const double threshold);

void bml_submatrix2matrix_ellpack_double_complex(
    const bml_matrix_dense_t * A,
    bml_matrix_ellpack_t * B,
    const int *core_halo_index,
    const int lsize,
    const int llsize,
    const double threshold);

void *bml_getVector_ellpack(
    const bml_matrix_ellpack_t * A,
    const int *jj,
    const int irow,
    const int colCnt);

void *bml_getVector_ellpack_single_real(
    const bml_matrix_ellpack_t * A,
    const int *jj,
    const int irow,
    const int colCnt);

void *bml_getVector_ellpack_double_real(
    const bml_matrix_ellpack_t * A,
    const int *jj,
    const int irow,
    const int colCnt);

void *bml_getVector_ellpack_single_complex(
    const bml_matrix_ellpack_t * A,
    const int *jj,
    const int irow,
    const int colCnt);

void *bml_getVector_ellpack_double_complex(
    const bml_matrix_ellpack_t * A,
    const int *jj,
    const int irow,
    const int colCnt);

bml_matrix_ellpack_t *bml_group_matrix_ellpack(
    const bml_matrix_ellpack_t * A,
    const int *hindex,
    const int ngroups,
    const double threshold);

bml_matrix_ellpack_t *bml_group_matrix_ellpack_single_real(
    const bml_matrix_ellpack_t * A,
    const int *hindex,
    const int ngroups,
    const double threshold);

bml_matrix_ellpack_t *bml_group_matrix_ellpack_double_real(
    const bml_matrix_ellpack_t * A,
    const int *hindex,
    const int ngroups,
    const double threshold);

bml_matrix_ellpack_t *bml_group_matrix_ellpack_single_complex(
    const bml_matrix_ellpack_t * A,
    const int *hindex,
    const int ngroups,
    const double threshold);

bml_matrix_ellpack_t *bml_group_matrix_ellpack_double_complex(
    const bml_matrix_ellpack_t * A,
    const int *hindex,
    const int ngroups,
    const double threshold);

void bml_adjacency_ellpack(
    const bml_matrix_ellpack_t * A,
    int *xadj,
    int *adjncy,
    const int base_flag);

void bml_adjacency_group_ellpack(
    const bml_matrix_ellpack_t * A,
    const int *hindex,
    const int nnodes,
    int *xadj,
    int *adjncy,
    const int base_flag);

#endif
