/** \file */

#ifndef __BML_ADD_H
#define __BML_ADD_H

#include "bml_types.h"

void bml_add(
    bml_matrix_t * A,
    bml_matrix_t * B,
    double alpha,
    double beta,
    double threshold);

double bml_add_norm(
    bml_matrix_t * A,
    bml_matrix_t * B,
    double alpha,
    double beta,
    double threshold);

void bml_add_identity(
    bml_matrix_t * A,
    double beta,
    double threshold);

void bml_scale_add_identity(
    bml_matrix_t * A,
    double alpha,
    double beta,
    double threshold);

#endif
