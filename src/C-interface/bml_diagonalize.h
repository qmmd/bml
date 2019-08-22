#ifndef __BML_DIAGONALIZE_H
#define __BML_DIAGONALIZE_H

#include "bml_types.h"

void bml_diagonalize(
    bml_matrix_t * A,
    void *eigenvalues,
    bml_matrix_t * eigenvectors);

#endif
