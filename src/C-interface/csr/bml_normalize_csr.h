#ifndef __BML_NORMALIZE_CSR_H
#define __BML_NORMALIZE_CSR_H

#include "bml_types_csr.h"

void bml_normalize_csr(
    bml_matrix_csr_t * A,
    double mineval,
    double maxeval);

void bml_normalize_csr_single_real(
    bml_matrix_csr_t * A,
    double mineval,
    double maxeval);

void bml_normalize_csr_double_real(
    bml_matrix_csr_t * A,
    double mineval,
    double maxeval);

void bml_normalize_csr_single_complex(
    bml_matrix_csr_t * A,
    double mineval,
    double maxeval);

void bml_normalize_csr_double_complex(
    bml_matrix_csr_t * A,
    double mineval,
    double maxeval);

void *bml_gershgorin_csr(
    bml_matrix_csr_t * A);

void *bml_gershgorin_csr_single_real(
    bml_matrix_csr_t * A);

void *bml_gershgorin_csr_double_real(
    bml_matrix_csr_t * A);

void *bml_gershgorin_csr_single_complex(
    bml_matrix_csr_t * A);

void *bml_gershgorin_csr_double_complex(
    bml_matrix_csr_t * A);

void *bml_gershgorin_partial_csr(
    bml_matrix_csr_t * A,
    int nrows);

void *bml_gershgorin_partial_csr_single_real(
    bml_matrix_csr_t * A,
    int nrows);

void *bml_gershgorin_partial_csr_double_real(
    bml_matrix_csr_t * A,
    int nrows);

void *bml_gershgorin_partial_csr_single_complex(
    bml_matrix_csr_t * A,
    int nrows);

void *bml_gershgorin_partial_csr_double_complex(
    bml_matrix_csr_t * A,
    int nrows);

#endif
