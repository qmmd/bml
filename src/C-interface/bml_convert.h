/** \file */

#ifndef __BML_CONVERT_H
#define __BML_CONVERT_H

#include "bml_types.h"

bml_matrix_t *bml_convert(
    bml_matrix_t * A,
    bml_matrix_type_t matrix_type,
    bml_matrix_precision_t matrix_precision,
    int M,
    bml_distribution_mode_t distrib_mode);

#endif
