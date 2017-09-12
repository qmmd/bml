#ifndef __TEST_BML_ADD_H
#define __TEST_BML_ADD_H

#include <bml.h>

int test_bml_add(
    const int N,
    const bml_matrix_type_t matrix_type,
    const bml_matrix_precision_t matrix_precision,
    const int M);

int test_bml_add_single_real(
    const int N,
    const bml_matrix_type_t matrix_type,
    const bml_matrix_precision_t matrix_precision,
    const int M);

int test_bml_add_double_real(
    const int N,
    const bml_matrix_type_t matrix_type,
    const bml_matrix_precision_t matrix_precision,
    const int M);

int test_bml_add_single_complex(
    const int N,
    const bml_matrix_type_t matrix_type,
    const bml_matrix_precision_t matrix_precision,
    const int M);

int test_bml_add_double_complex(
    const int N,
    const bml_matrix_type_t matrix_type,
    const bml_matrix_precision_t matrix_precision,
    const int M);

#endif
