#ifndef __BML_THRESHOLD_CSR_H
#define __BML_THRESHOLD_CSR_H

#include "bml_types_csr.h"

bml_matrix_csr_t *bml_threshold_new_csr(
    bml_matrix_csr_t * A,
    double threshold);

bml_matrix_csr_t
    * bml_threshold_new_csr_single_real(bml_matrix_csr_t * A,
                                        double threshold);

bml_matrix_csr_t
    * bml_threshold_new_csr_double_real(bml_matrix_csr_t * A,
                                        double threshold);

bml_matrix_csr_t
    * bml_threshold_new_csr_single_complex(bml_matrix_csr_t * A,
                                           double threshold);

bml_matrix_csr_t
    * bml_threshold_new_csr_double_complex(bml_matrix_csr_t * A,
                                           double threshold);

void bml_threshold_csr(
    bml_matrix_csr_t * A,
    double threshold);

void bml_threshold_csr_single_real(
    bml_matrix_csr_t * A,
    double threshold);

void bml_threshold_csr_double_real(
    bml_matrix_csr_t * A,
    double threshold);

void bml_threshold_csr_single_complex(
    bml_matrix_csr_t * A,
    double threshold);

void bml_threshold_csr_double_complex(
    bml_matrix_csr_t * A,
    double threshold);

#endif
