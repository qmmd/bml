#ifndef __BML_INVERSE_ELLPACK_H
#define __BML_INVERSE_ELLPACK_H

#include "bml_types_ellpack.h"

bml_matrix_ellpack_t *bml_inverse_ellpack(
    bml_matrix_ellpack_t * A);

bml_matrix_ellpack_t *bml_inverse_ellpack_single_real(
    bml_matrix_ellpack_t * A);

bml_matrix_ellpack_t *bml_inverse_ellpack_double_real(
    bml_matrix_ellpack_t * A);

bml_matrix_ellpack_t *bml_inverse_ellpack_single_complex(
    bml_matrix_ellpack_t * A);

bml_matrix_ellpack_t *bml_inverse_ellpack_double_complex(
    bml_matrix_ellpack_t * A);

#endif
