#ifndef __BML_INVERSE_CSR_H
#define __BML_INVERSE_CSR_H

#include "bml_types_csr.h"

bml_matrix_csr_t *bml_inverse_csr(
    bml_matrix_csr_t * A);

bml_matrix_csr_t *bml_inverse_csr_single_real(
    bml_matrix_csr_t * A);

bml_matrix_csr_t *bml_inverse_csr_double_real(
    bml_matrix_csr_t * A);

bml_matrix_csr_t *bml_inverse_csr_single_complex(
    bml_matrix_csr_t * A);

bml_matrix_csr_t *bml_inverse_csr_double_complex(
    bml_matrix_csr_t * A);

#endif
