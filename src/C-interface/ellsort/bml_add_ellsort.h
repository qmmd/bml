#ifndef __BML_ADD_ELLSORT_H
#define __BML_ADD_ELLSORT_H

#include "bml_types_ellsort.h"

void bml_add_ellsort(
    bml_matrix_ellsort_t * A,
    bml_matrix_ellsort_t * B,
    double alpha,
    double beta,
    double threshold);

void bml_add_ellsort_single_real(
    bml_matrix_ellsort_t * A,
    bml_matrix_ellsort_t * B,
    double alpha,
    double beta,
    double threshold);

void bml_add_ellsort_double_real(
    bml_matrix_ellsort_t * A,
    bml_matrix_ellsort_t * B,
    double alpha,
    double beta,
    double threshold);

void bml_add_ellsort_single_complex(
    bml_matrix_ellsort_t * A,
    bml_matrix_ellsort_t * B,
    double alpha,
    double beta,
    double threshold);

void bml_add_ellsort_double_complex(
    bml_matrix_ellsort_t * A,
    bml_matrix_ellsort_t * B,
    double alpha,
    double beta,
    double threshold);

double bml_add_norm_ellsort(
    bml_matrix_ellsort_t * A,
    bml_matrix_ellsort_t * B,
    double alpha,
    double beta,
    double threshold);

double bml_add_norm_ellsort_single_real(
    bml_matrix_ellsort_t * A,
    bml_matrix_ellsort_t * B,
    double alpha,
    double beta,
    double threshold);

double bml_add_norm_ellsort_double_real(
    bml_matrix_ellsort_t * A,
    bml_matrix_ellsort_t * B,
    double alpha,
    double beta,
    double threshold);

double bml_add_norm_ellsort_single_complex(
    bml_matrix_ellsort_t * A,
    bml_matrix_ellsort_t * B,
    double alpha,
    double beta,
    double threshold);

double bml_add_norm_ellsort_double_complex(
    bml_matrix_ellsort_t * A,
    bml_matrix_ellsort_t * B,
    double alpha,
    double beta,
    double threshold);

void bml_add_identity_ellsort(
    bml_matrix_ellsort_t * A,
    double beta,
    double threshold);

void bml_add_identity_ellsort_single_real(
    bml_matrix_ellsort_t * A,
    double beta,
    double threshold);

void bml_add_identity_ellsort_double_real(
    bml_matrix_ellsort_t * A,
    double beta,
    double threshold);

void bml_add_identity_ellsort_single_complex(
    bml_matrix_ellsort_t * A,
    double beta,
    double threshold);

void bml_add_identity_ellsort_double_complex(
    bml_matrix_ellsort_t * A,
    double beta,
    double threshold);

void bml_scale_add_identity_ellsort(
    bml_matrix_ellsort_t * A,
    double alpha,
    double beta,
    double threshold);

void bml_scale_add_identity_ellsort_single_real(
    bml_matrix_ellsort_t * A,
    double alpha,
    double beta,
    double threshold);

void bml_scale_add_identity_ellsort_double_real(
    bml_matrix_ellsort_t * A,
    double alpha,
    double beta,
    double threshold);

void bml_scale_add_identity_ellsort_single_complex(
    bml_matrix_ellsort_t * A,
    double alpha,
    double beta,
    double threshold);

void bml_scale_add_identity_ellsort_double_complex(
    bml_matrix_ellsort_t * A,
    double alpha,
    double beta,
    double threshold);

#endif
