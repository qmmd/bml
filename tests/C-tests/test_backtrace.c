#include "bml_logger.h"

void
foo(
    void)
{
    LOG_DEBUG("Test function\n");
    LOG_ERROR("Print backtrace\n");
}

int
main(
    int argc,
    char **argv)
{
    LOG_INFO("Starting test...\n");
    foo();
}
