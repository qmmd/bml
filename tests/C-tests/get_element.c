#include "bml.h"
#include "bml_test.h"

#include <stdio.h>

int
test_get_element(
    const int N,
    const bml_matrix_type_t matrix_type,
    const bml_matrix_precision_t matrix_precision,
    const int M)
{
    switch (matrix_precision)
    {
        case single_real:
            return test_get_element_single_real(N, matrix_type,
                                                matrix_precision, M);
            break;
        case double_real:
            return test_get_element_double_real(N, matrix_type,
                                                matrix_precision, M);
            break;
        case single_complex:
            return test_get_element_single_complex(N, matrix_type,
                                                   matrix_precision, M);
            break;
        case double_complex:
            return test_get_element_double_complex(N, matrix_type,
                                                   matrix_precision, M);
            break;
        default:
            fprintf(stderr, "unknown matrix precision\n");
            return -1;
            break;
    }
}
