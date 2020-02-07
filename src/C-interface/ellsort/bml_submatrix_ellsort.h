#ifndef __BML_SUBMATRIX_ELLSORT_H
#define __BML_SUBMATRIX_ELLSORT_H

#include "bml_types_ellsort.h"
#include "../dense/bml_types_dense.h"

#include <complex.h>

void bml_matrix2submatrix_index_ellsort(
    bml_matrix_ellsort_t * A,
    bml_matrix_ellsort_t * B,
    int *nodelist,
    int nsize,
    int *core_halo_index,
    int *vsize,
    int double_jump_flag);

void bml_matrix2submatrix_index_ellsort_single_real(
    bml_matrix_ellsort_t * A,
    bml_matrix_ellsort_t * B,
    int *nodelist,
    int nsize,
    int *core_halo_index,
    int *vsize,
    int double_jump_flag);

void bml_matrix2submatrix_index_ellsort_double_real(
    bml_matrix_ellsort_t * A,
    bml_matrix_ellsort_t * B,
    int *nodelist,
    int nsize,
    int *core_halo_index,
    int *vsize,
    int double_jump_flag);

void bml_matrix2submatrix_index_ellsort_single_complex(
    bml_matrix_ellsort_t * A,
    bml_matrix_ellsort_t * B,
    int *nodelist,
    int nsize,
    int *core_halo_index,
    int *vsize,
    int double_jump_flag);

void bml_matrix2submatrix_index_ellsort_double_complex(
    bml_matrix_ellsort_t * A,
    bml_matrix_ellsort_t * B,
    int *nodelist,
    int nsize,
    int *core_halo_index,
    int *vsize,
    int double_jump_flag);

void bml_matrix2submatrix_index_graph_ellsort(
    bml_matrix_ellsort_t * B,
    int *nodelist,
    int nsize,
    int *core_halo_index,
    int *vsize,
    int double_jump_flag);

void bml_matrix2submatrix_index_graph_ellsort_single_real(
    bml_matrix_ellsort_t * B,
    int *nodelist,
    int nsize,
    int *core_halo_index,
    int *vsize,
    int double_jump_flag);

void bml_matrix2submatrix_index_graph_ellsort_double_real(
    bml_matrix_ellsort_t * B,
    int *nodelist,
    int nsize,
    int *core_halo_index,
    int *vsize,
    int double_jump_flag);

void bml_matrix2submatrix_index_graph_ellsort_single_complex(
    bml_matrix_ellsort_t * B,
    int *nodelist,
    int nsize,
    int *core_halo_index,
    int *vsize,
    int double_jump_flag);

void bml_matrix2submatrix_index_graph_ellsort_double_complex(
    bml_matrix_ellsort_t * B,
    int *nodelist,
    int nsize,
    int *core_halo_index,
    int *vsize,
    int double_jump_flag);

void bml_matrix2submatrix_ellsort(
    bml_matrix_ellsort_t * A,
    bml_matrix_dense_t * B,
    int *core_halo_index,
    int lsize);

void bml_matrix2submatrix_ellsort_single_real(
    bml_matrix_ellsort_t * A,
    bml_matrix_dense_t * B,
    int *core_halo_index,
    int lsize);

void bml_matrix2submatrix_ellsort_double_real(
    bml_matrix_ellsort_t * A,
    bml_matrix_dense_t * B,
    int *core_halo_index,
    int lsize);

void bml_matrix2submatrix_ellsort_single_complex(
    bml_matrix_ellsort_t * A,
    bml_matrix_dense_t * B,
    int *core_halo_index,
    int lsize);

void bml_matrix2submatrix_ellsort_double_complex(
    bml_matrix_ellsort_t * A,
    bml_matrix_dense_t * B,
    int *core_halo_index,
    int lsize);

void bml_submatrix2matrix_ellsort(
    bml_matrix_dense_t * A,
    bml_matrix_ellsort_t * B,
    int *core_halo_index,
    int lsize,
    int llsize,
    double threshold);

void bml_submatrix2matrix_ellsort_single_real(
    bml_matrix_dense_t * A,
    bml_matrix_ellsort_t * B,
    int *core_halo_index,
    int lsize,
    int llsize,
    double threshold);

void bml_submatrix2matrix_ellsort_double_real(
    bml_matrix_dense_t * A,
    bml_matrix_ellsort_t * B,
    int *core_halo_index,
    int lsize,
    int llsize,
    double threshold);

void bml_submatrix2matrix_ellsort_single_complex(
    bml_matrix_dense_t * A,
    bml_matrix_ellsort_t * B,
    int *core_halo_index,
    int lsize,
    int llsize,
    double threshold);

void bml_submatrix2matrix_ellsort_double_complex(
    bml_matrix_dense_t * A,
    bml_matrix_ellsort_t * B,
    int *core_halo_index,
    int lsize,
    int llsize,
    double threshold);

void *bml_getVector_ellsort(
    bml_matrix_ellsort_t * A,
    int *jj,
    int irow,
    int colCnt);

void *bml_getVector_ellsort_single_real(
    bml_matrix_ellsort_t * A,
    int *jj,
    int irow,
    int colCnt);

void *bml_getVector_ellsort_double_real(
    bml_matrix_ellsort_t * A,
    int *jj,
    int irow,
    int colCnt);

void *bml_getVector_ellsort_single_complex(
    bml_matrix_ellsort_t * A,
    int *jj,
    int irow,
    int colCnt);

void *bml_getVector_ellsort_double_complex(
    bml_matrix_ellsort_t * A,
    int *jj,
    int irow,
    int colCnt);

bml_matrix_ellsort_t *bml_group_matrix_ellsort(
    bml_matrix_ellsort_t * A,
    int *hindex,
    int ngroups,
    double threshold);

bml_matrix_ellsort_t
    * bml_group_matrix_ellsort_single_real(bml_matrix_ellsort_t * A,
                                           int *hindex,
                                           int ngroups,
                                           double threshold);

bml_matrix_ellsort_t
    * bml_group_matrix_ellsort_double_real(bml_matrix_ellsort_t * A,
                                           int *hindex,
                                           int ngroups,
                                           double threshold);

bml_matrix_ellsort_t
    * bml_group_matrix_ellsort_single_complex(bml_matrix_ellsort_t * A,
                                              int *hindex,
                                              int ngroups,
                                              double threshold);

bml_matrix_ellsort_t
    * bml_group_matrix_ellsort_double_complex(bml_matrix_ellsort_t * A,
                                              int *hindex,
                                              int ngroups,
                                              double threshold);

void bml_adjacency_ellsort(
    bml_matrix_ellsort_t * A,
    int *xadj,
    int *adjncy,
    int base_flag);

void bml_adjacency_group_ellsort(
    bml_matrix_ellsort_t * A,
    int *hindex,
    int nnodes,
    int *xadj,
    int *adjncy,
    int base_flag);

#endif
