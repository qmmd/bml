/** \file */

#ifndef __BML_CONVERT_H
#define __BML_CONVERT_H

#include "bml_types.h"

bml_matrix_t *bml_convert(
    const bml_matrix_t * A,
    const bml_matrix_type_t matrix_type,
    const bml_matrix_precision_t matrix_precision,
    const int M,
    const bml_distribution_mode_t distrib_mode);

#endif
