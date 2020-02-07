#ifndef __BML_THRESHOLD_ELLPACK_H
#define __BML_THRESHOLD_ELLPACK_H

#include "bml_types_ellpack.h"

bml_matrix_ellpack_t *bml_threshold_new_ellpack(
    bml_matrix_ellpack_t * A,
    double threshold);

bml_matrix_ellpack_t
    * bml_threshold_new_ellpack_single_real(bml_matrix_ellpack_t * A,
                                            double threshold);

bml_matrix_ellpack_t
    * bml_threshold_new_ellpack_double_real(bml_matrix_ellpack_t * A,
                                            double threshold);

bml_matrix_ellpack_t
    * bml_threshold_new_ellpack_single_complex(bml_matrix_ellpack_t * A,
                                               double threshold);

bml_matrix_ellpack_t
    * bml_threshold_new_ellpack_double_complex(bml_matrix_ellpack_t * A,
                                               double threshold);

void bml_threshold_ellpack(
    bml_matrix_ellpack_t * A,
    double threshold);

void bml_threshold_ellpack_single_real(
    bml_matrix_ellpack_t * A,
    double threshold);

void bml_threshold_ellpack_double_real(
    bml_matrix_ellpack_t * A,
    double threshold);

void bml_threshold_ellpack_single_complex(
    bml_matrix_ellpack_t * A,
    double threshold);

void bml_threshold_ellpack_double_complex(
    bml_matrix_ellpack_t * A,
    double threshold);

#endif
