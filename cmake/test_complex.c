#include <complex.h>

int
main(
    int argc,
    char **argv)
{
    double complex x;

    x = 3.0 + 4.0 * _Complex_I;
    if (creal(x) == 3.0)
    {
        return 0;
    }
    return 1;
}
