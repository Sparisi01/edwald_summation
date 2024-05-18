#ifndef COMPLEX_UTILS_H
#define COMPLEX_UTILS_H

#include <complex.h>
#include <stdio.h>
#include <tgmath.h>

void cprint(double complex z)
{
    printf("%.1f%+.1fi", creal(z), cimag(z));
}

void cprintln(double complex z)
{
    cprint(z);
    printf("\n");
}

#endif