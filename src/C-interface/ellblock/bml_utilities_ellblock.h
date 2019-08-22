/** \file */

#ifndef __BML_UTILITIES_ELLBLOCK_H
#define __BML_UTILITIES_ELLBLOCK_H

#include "bml_types_ellblock.h"

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

void bml_read_bml_matrix_ellblock(
    bml_matrix_ellblock_t * A,
    char *filename);

void bml_read_bml_matrix_ellblock_single_real(
    bml_matrix_ellblock_t * A,
    char *filename);

void bml_read_bml_matrix_ellblock_double_real(
    bml_matrix_ellblock_t * A,
    char *filename);

void bml_read_bml_matrix_ellblock_single_complex(
    bml_matrix_ellblock_t * A,
    char *filename);

void bml_read_bml_matrix_ellblock_double_complex(
    bml_matrix_ellblock_t * A,
    char *filename);

void bml_write_bml_matrix_ellblock(
    bml_matrix_ellblock_t * A,
    char *filename);

void bml_write_bml_matrix_ellblock_single_real(
    bml_matrix_ellblock_t * A,
    char *filename);

void bml_write_bml_matrix_ellblock_double_real(
    bml_matrix_ellblock_t * A,
    char *filename);

void bml_write_bml_matrix_ellblock_single_complex(
    bml_matrix_ellblock_t * A,
    char *filename);

void bml_write_bml_matrix_ellblock_double_complex(
    bml_matrix_ellblock_t * A,
    char *filename);

double bml_sum_squares_single_real(
    void *v,
    int,
    int,
    int);
double bml_sum_squares_double_real(
    void *v,
    int,
    int,
    int);
double bml_sum_squares_single_complex(
    void *v,
    int,
    int,
    int);
double bml_sum_squares_double_complex(
    void *v,
    int,
    int,
    int);

double bml_norm_inf_single_real(
    void *v,
    int,
    int,
    int);
double bml_norm_inf_double_real(
    void *v,
    int,
    int,
    int);
double bml_norm_inf_single_complex(
    void *v,
    int,
    int,
    int);
double bml_norm_inf_double_complex(
    void *v,
    int,
    int,
    int);

double bml_norm_inf_fast_single_real(
    void *v,
    int);
double bml_norm_inf_fast_double_real(
    void *v,
    int);
double bml_norm_inf_fast_single_complex(
    void *v,
    int);
double bml_norm_inf_fast_double_complex(
    void *v,
    int);

#endif
