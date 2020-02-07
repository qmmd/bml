#ifndef __BML_TRANSPOSE_ELLPACK_H
#define __BML_TRANSPOSE_ELLPACK_H

#include "bml_types_ellpack.h"

bml_matrix_ellpack_t *bml_transpose_new_ellpack(
    bml_matrix_ellpack_t * A);

bml_matrix_ellpack_t
    * bml_transpose_new_ellpack_single_real(bml_matrix_ellpack_t * A);

bml_matrix_ellpack_t
    * bml_transpose_new_ellpack_double_real(bml_matrix_ellpack_t * A);

bml_matrix_ellpack_t
    * bml_transpose_new_ellpack_single_complex(bml_matrix_ellpack_t * A);

bml_matrix_ellpack_t
    * bml_transpose_new_ellpack_double_complex(bml_matrix_ellpack_t * A);

void bml_transpose_ellpack(
    bml_matrix_ellpack_t * A);

void bml_transpose_ellpack_single_real(
    bml_matrix_ellpack_t * A);

void bml_transpose_ellpack_double_real(
    bml_matrix_ellpack_t * A);

void bml_transpose_ellpack_single_complex(
    bml_matrix_ellpack_t * A);

void bml_transpose_ellpack_double_complex(
    bml_matrix_ellpack_t * A);

#endif
