#include "bml_init.h"
#include "bml_parallel.h"

#ifdef BML_USE_MAGMA
#include "magma_v2.h"
#endif
#ifdef BML_USE_XSMM
#include "libxsmm.h"
#endif

#include <stdlib.h>

/** Initialize.
 *
 * \ingroup init_group_C
 *
 * \param argc Number of args
 * \param argv Args
 */
void
bml_init(
#ifdef DO_MPI
    MPI_Comm comm
#endif
    )
{
#ifdef BML_USE_MAGMA
    magma_init();
#endif
#ifdef BML_USE_XSMM
    // Initialize the library; pay for setup cost here
    libxsmm_init();
#endif
#ifdef DO_MPI
    bml_initParallel(comm);
#endif
}

/** Initialize from Fortran.
 *
 * \ingroup init_group_C
 *
 * \param Comm from Fortran
 */
void
bml_initF(
    int fcomm)
{
    bml_initParallelF(fcomm);
}
