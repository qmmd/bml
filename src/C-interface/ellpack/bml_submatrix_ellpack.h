#ifndef __BML_SUBMATRIX_ELLPACK_H
#define __BML_SUBMATRIX_ELLPACK_H

#include "bml_types_ellpack.h"
#include "../dense/bml_types_dense.h"

#include <complex.h>

void bml_matrix2submatrix_index_ellpack(
    bml_matrix_ellpack_t * A,
    bml_matrix_ellpack_t * B,
    int *nodelist,
    int nsize,
    int *core_halo_index,
    int *vsize,
    int double_jump_flag);

void bml_matrix2submatrix_index_ellpack_single_real(
    bml_matrix_ellpack_t * A,
    bml_matrix_ellpack_t * B,
    int *nodelist,
    int nsize,
    int *core_halo_index,
    int *vsize,
    int double_jump_flag);

void bml_matrix2submatrix_index_ellpack_double_real(
    bml_matrix_ellpack_t * A,
    bml_matrix_ellpack_t * B,
    int *nodelist,
    int nsize,
    int *core_halo_index,
    int *vsize,
    int double_jump_flag);

void bml_matrix2submatrix_index_ellpack_single_complex(
    bml_matrix_ellpack_t * A,
    bml_matrix_ellpack_t * B,
    int *nodelist,
    int nsize,
    int *core_halo_index,
    int *vsize,
    int double_jump_flag);

void bml_matrix2submatrix_index_ellpack_double_complex(
    bml_matrix_ellpack_t * A,
    bml_matrix_ellpack_t * B,
    int *nodelist,
    int nsize,
    int *core_halo_index,
    int *vsize,
    int double_jump_flag);

void bml_matrix2submatrix_index_graph_ellpack(
    bml_matrix_ellpack_t * B,
    int *nodelist,
    int nsize,
    int *core_halo_index,
    int *vsize,
    int double_jump_flag);

void bml_matrix2submatrix_index_graph_ellpack_single_real(
    bml_matrix_ellpack_t * B,
    int *nodelist,
    int nsize,
    int *core_halo_index,
    int *vsize,
    int double_jump_flag);

void bml_matrix2submatrix_index_graph_ellpack_double_real(
    bml_matrix_ellpack_t * B,
    int *nodelist,
    int nsize,
    int *core_halo_index,
    int *vsize,
    int double_jump_flag);

void bml_matrix2submatrix_index_graph_ellpack_single_complex(
    bml_matrix_ellpack_t * B,
    int *nodelist,
    int nsize,
    int *core_halo_index,
    int *vsize,
    int double_jump_flag);

void bml_matrix2submatrix_index_graph_ellpack_double_complex(
    bml_matrix_ellpack_t * B,
    int *nodelist,
    int nsize,
    int *core_halo_index,
    int *vsize,
    int double_jump_flag);

void bml_matrix2submatrix_ellpack(
    bml_matrix_ellpack_t * A,
    bml_matrix_dense_t * B,
    int *core_halo_index,
    int lsize);

void bml_matrix2submatrix_ellpack_single_real(
    bml_matrix_ellpack_t * A,
    bml_matrix_dense_t * B,
    int *core_halo_index,
    int lsize);

void bml_matrix2submatrix_ellpack_double_real(
    bml_matrix_ellpack_t * A,
    bml_matrix_dense_t * B,
    int *core_halo_index,
    int lsize);

void bml_matrix2submatrix_ellpack_single_complex(
    bml_matrix_ellpack_t * A,
    bml_matrix_dense_t * B,
    int *core_halo_index,
    int lsize);

void bml_matrix2submatrix_ellpack_double_complex(
    bml_matrix_ellpack_t * A,
    bml_matrix_dense_t * B,
    int *core_halo_index,
    int lsize);

void bml_submatrix2matrix_ellpack(
    bml_matrix_dense_t * A,
    bml_matrix_ellpack_t * B,
    int *core_halo_index,
    int lsize,
    int llsize,
    double threshold);

void bml_submatrix2matrix_ellpack_single_real(
    bml_matrix_dense_t * A,
    bml_matrix_ellpack_t * B,
    int *core_halo_index,
    int lsize,
    int llsize,
    double threshold);

void bml_submatrix2matrix_ellpack_double_real(
    bml_matrix_dense_t * A,
    bml_matrix_ellpack_t * B,
    int *core_halo_index,
    int lsize,
    int llsize,
    double threshold);

void bml_submatrix2matrix_ellpack_single_complex(
    bml_matrix_dense_t * A,
    bml_matrix_ellpack_t * B,
    int *core_halo_index,
    int lsize,
    int llsize,
    double threshold);

void bml_submatrix2matrix_ellpack_double_complex(
    bml_matrix_dense_t * A,
    bml_matrix_ellpack_t * B,
    int *core_halo_index,
    int lsize,
    int llsize,
    double threshold);

void *bml_getVector_ellpack(
    bml_matrix_ellpack_t * A,
    int *jj,
    int irow,
    int colCnt);

void *bml_getVector_ellpack_single_real(
    bml_matrix_ellpack_t * A,
    int *jj,
    int irow,
    int colCnt);

void *bml_getVector_ellpack_double_real(
    bml_matrix_ellpack_t * A,
    int *jj,
    int irow,
    int colCnt);

void *bml_getVector_ellpack_single_complex(
    bml_matrix_ellpack_t * A,
    int *jj,
    int irow,
    int colCnt);

void *bml_getVector_ellpack_double_complex(
    bml_matrix_ellpack_t * A,
    int *jj,
    int irow,
    int colCnt);

bml_matrix_ellpack_t *bml_group_matrix_ellpack(
    bml_matrix_ellpack_t * A,
    int *hindex,
    int ngroups,
    double threshold);

bml_matrix_ellpack_t
    * bml_group_matrix_ellpack_single_real(bml_matrix_ellpack_t * A,
                                           int *hindex,
                                           int ngroups,
                                           double threshold);

bml_matrix_ellpack_t
    * bml_group_matrix_ellpack_double_real(bml_matrix_ellpack_t * A,
                                           int *hindex,
                                           int ngroups,
                                           double threshold);

bml_matrix_ellpack_t
    * bml_group_matrix_ellpack_single_complex(bml_matrix_ellpack_t * A,
                                              int *hindex,
                                              int ngroups,
                                              double threshold);

bml_matrix_ellpack_t
    * bml_group_matrix_ellpack_double_complex(bml_matrix_ellpack_t * A,
                                              int *hindex,
                                              int ngroups,
                                              double threshold);

void bml_adjacency_ellpack(
    bml_matrix_ellpack_t * A,
    int *xadj,
    int *adjncy,
    int base_flag);

void bml_adjacency_group_ellpack(
    bml_matrix_ellpack_t * A,
    int *hindex,
    int nnodes,
    int *xadj,
    int *adjncy,
    int base_flag);

#endif
