#ifndef __BML_TEST_H
#define __BML_TEST_H

#include <bml.h>

typedef int (
    *test_function_t) (
    const int N,
    const bml_matrix_type_t matrix_type,
    const bml_matrix_precision_t matrix_precision,
    const int M);

#include "add_matrix.h"
#include "adjacency_matrix.h"
#include "adjungate_triangle_matrix.h"
#include "allocate_matrix.h"
#include "convert_matrix.h"
#include "copy_matrix.h"
#include "diagonalize_matrix.h"
#include "get_element.h"
#include "get_set_diagonal.h"
#include "get_sparsity.h"
#include "import_export_matrix.h"
#include "inverse_matrix.h"
#include "multiply_matrix.h"
#include "normalize_matrix.h"
#include "norm_matrix.h"
#include "scale_matrix.h"
#include "set_row.h"
#include "submatrix_matrix.h"
#include "test_bml_gemm.h"
#include "threshold_matrix.h"
#include "trace_matrix.h"
#include "transpose_matrix.h"

#endif
