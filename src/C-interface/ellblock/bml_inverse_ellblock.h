#ifndef __BML_INVERSE_ELLBLOCK_H
#define __BML_INVERSE_ELLBLOCK_H

#include "bml_types_ellblock.h"

bml_matrix_ellblock_t *bml_inverse_ellblock(
    bml_matrix_ellblock_t * A);

bml_matrix_ellblock_t *bml_inverse_ellblock_single_real(
    bml_matrix_ellblock_t * A);

bml_matrix_ellblock_t *bml_inverse_ellblock_double_real(
    bml_matrix_ellblock_t * A);

bml_matrix_ellblock_t *bml_inverse_ellblock_single_complex(
    bml_matrix_ellblock_t * A);

bml_matrix_ellblock_t *bml_inverse_ellblock_double_complex(
    bml_matrix_ellblock_t * A);

#endif
