#include "bml.h"
#include "bml_test.h"

#define STRINGIFY2(a) #a
#define STRINGIFY(a) STRINGIFY2(a)

int
main(
    int argc,
    char **argv)
{
    const int N = 7;
    const int M = 7;

    bml_log(BML_LOG_INFO, "testing %s:%s\n",
            STRINGIFY(MATRIX_TYPE_NAME), STRINGIFY(MATRIX_PRECISION));
    if (test_function(N, MATRIX_TYPE_NAME, MATRIX_PRECISION, M) != 0)
    {
        LOG_ERROR("test failed\n");
        return -1;
    }
    LOG_INFO("test passed\n");
    return 0;
}
