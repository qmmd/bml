#include "bml_shutdown.h"
#include "bml_parallel.h"

#ifdef BML_USE_MAGMA
#include "magma_v2.h"
#endif
#ifdef BML_USE_XSMM
#include "libxsmm.h"
#endif

#include <stdlib.h>

/** Shutdown.
 *
 * \ingroup shutdown_group_C
 *
 */
void
bml_shutdown(
    )
{
    bml_shutdownParallel();
#ifdef BML_USE_XSMM
    // De-initialize the library and free internal memory (optional)
    libxsmm_finalize();
#endif
#ifdef BML_USE_MAGMA
    magma_finalize();
#endif
}

/** Shutdown from Fortran.
 *
 * \ingroup shutdown_group_C
 *
 */
void
bml_shutdownF(
    )
{
    bml_shutdownParallelF();
    // Future: shutdown GPUs, cublas, cusparse, etc.
}
