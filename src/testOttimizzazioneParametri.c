#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <tgmath.h>

#include "includes/coulomb/coulomb_potential.h"
#include "includes/ewald/edwald.h"
#include "includes/utils/statistic.h"

int _N_PARTICLES = 4000;
double _CELL_LENGHT = 1;

int main(int argc, char const *argv[])
{
    double Q = 0.01 * _N_PARTICLES;

    optimizeParameter(1e-4, _CELL_LENGHT, _N_PARTICLES, Q);

    return 0;
}
