#ifndef __BML_NORMALIZE_BELLBLOCK_H
#define __BML_NORMALIZE_BELLBLOCK_H

#include "bml_types_ellblock.h"

void bml_normalize_ellblock(
    bml_matrix_ellblock_t * A,
    double mineval,
    double maxeval);

void bml_normalize_ellblock_single_real(
    bml_matrix_ellblock_t * A,
    double mineval,
    double maxeval);

void bml_normalize_ellblock_double_real(
    bml_matrix_ellblock_t * A,
    double mineval,
    double maxeval);

void bml_normalize_ellblock_single_complex(
    bml_matrix_ellblock_t * A,
    double mineval,
    double maxeval);

void bml_normalize_ellblock_double_complex(
    bml_matrix_ellblock_t * A,
    double mineval,
    double maxeval);

void *bml_gershgorin_ellblock(
    bml_matrix_ellblock_t * A);

void *bml_gershgorin_ellblock_single_real(
    bml_matrix_ellblock_t * A);

void *bml_gershgorin_ellblock_double_real(
    bml_matrix_ellblock_t * A);

void *bml_gershgorin_ellblock_single_complex(
    bml_matrix_ellblock_t * A);

void *bml_gershgorin_ellblock_double_complex(
    bml_matrix_ellblock_t * A);

#endif
