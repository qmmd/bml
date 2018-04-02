#include "../ellsort/bml_adjungate_triangle_ellsort.h"
#include "bml_adjungate_triangle.h"
#include "bml_introspection.h"
#include "bml_logger.h"
#include "bml_types.h"
#include "bml_types_ellsort.h"

#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

void
bml_adjungate_triangle_ellsort(
    bml_matrix_ellsort_t * A,
    char *triangle)
{
//   LOG_ERROR("unkonwn precision for bml_adjungate_triangle\n");

    switch (bml_get_precision(A))
    {
        case single_real:
            bml_adjungate_triangle_ellsort_single_real(A, triangle);
            break;
        case double_real:
            bml_adjungate_triangle_ellsort_double_real(A, triangle);
            break;
        case single_complex:
            bml_adjungate_triangle_ellsort_single_complex(A, triangle);
            break;
        case double_complex:
            bml_adjungate_triangle_ellsort_double_complex(A, triangle);
            break;
        default:
            LOG_ERROR("unkonwn precision for bml_adjungate_triangle\n");
            break;
    }
}
