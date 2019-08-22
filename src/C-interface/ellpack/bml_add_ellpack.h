#ifndef __BML_ADD_ELLPACK_H
#define __BML_ADD_ELLPACK_H

#include "bml_types_ellpack.h"

void bml_add_ellpack(
    bml_matrix_ellpack_t * A,
    bml_matrix_ellpack_t * B,
    double alpha,
    double beta,
    double threshold);

void bml_add_ellpack_single_real(
    bml_matrix_ellpack_t * A,
    bml_matrix_ellpack_t * B,
    double alpha,
    double beta,
    double threshold);

void bml_add_ellpack_double_real(
    bml_matrix_ellpack_t * A,
    bml_matrix_ellpack_t * B,
    double alpha,
    double beta,
    double threshold);

void bml_add_ellpack_single_complex(
    bml_matrix_ellpack_t * A,
    bml_matrix_ellpack_t * B,
    double alpha,
    double beta,
    double threshold);

void bml_add_ellpack_double_complex(
    bml_matrix_ellpack_t * A,
    bml_matrix_ellpack_t * B,
    double alpha,
    double beta,
    double threshold);

double bml_add_norm_ellpack(
    bml_matrix_ellpack_t * A,
    bml_matrix_ellpack_t * B,
    double alpha,
    double beta,
    double threshold);

double bml_add_norm_ellpack_single_real(
    bml_matrix_ellpack_t * A,
    bml_matrix_ellpack_t * B,
    double alpha,
    double beta,
    double threshold);

double bml_add_norm_ellpack_double_real(
    bml_matrix_ellpack_t * A,
    bml_matrix_ellpack_t * B,
    double alpha,
    double beta,
    double threshold);

double bml_add_norm_ellpack_single_complex(
    bml_matrix_ellpack_t * A,
    bml_matrix_ellpack_t * B,
    double alpha,
    double beta,
    double threshold);

double bml_add_norm_ellpack_double_complex(
    bml_matrix_ellpack_t * A,
    bml_matrix_ellpack_t * B,
    double alpha,
    double beta,
    double threshold);

void bml_add_identity_ellpack(
    bml_matrix_ellpack_t * A,
    double beta,
    double threshold);

void bml_add_identity_ellpack_single_real(
    bml_matrix_ellpack_t * A,
    double beta,
    double threshold);

void bml_add_identity_ellpack_double_real(
    bml_matrix_ellpack_t * A,
    double beta,
    double threshold);

void bml_add_identity_ellpack_single_complex(
    bml_matrix_ellpack_t * A,
    double beta,
    double threshold);

void bml_add_identity_ellpack_double_complex(
    bml_matrix_ellpack_t * A,
    double beta,
    double threshold);

void bml_scale_add_identity_ellpack(
    bml_matrix_ellpack_t * A,
    double alpha,
    double beta,
    double threshold);

void bml_scale_add_identity_ellpack_single_real(
    bml_matrix_ellpack_t * A,
    double alpha,
    double beta,
    double threshold);

void bml_scale_add_identity_ellpack_double_real(
    bml_matrix_ellpack_t * A,
    double alpha,
    double beta,
    double threshold);

void bml_scale_add_identity_ellpack_single_complex(
    bml_matrix_ellpack_t * A,
    double alpha,
    double beta,
    double threshold);

void bml_scale_add_identity_ellpack_double_complex(
    bml_matrix_ellpack_t * A,
    double alpha,
    double beta,
    double threshold);

#endif
