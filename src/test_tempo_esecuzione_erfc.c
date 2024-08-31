#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <tgmath.h>
#include <time.h>

#include "includes/coulomb/coulomb_potential.h"
#include "includes/ewald/edwald.h"
#include "includes/utils/statistic.h"

int _N_PARTICLES = 100;
double _CELL_LENGHT = 1;
double _SIGMA_VELOCITIES = 1.;
double _CHARGE = 1.;

int main(void)
{

    clock_t c0, c1;

    c0 = clock();
    int N = 0;
    float dx = 10.0f - nextafterf(10, 0);
    for (float x = 10.0f; x >= 0.0f; x -= dx)
    {
        N++;
        erfc(x);
    }
    c1 = clock();
    printf("%d time: %g\n", N, (double)(c1 - c0) / CLOCKS_PER_SEC);

    c0 = clock();
    N = 0;
    dx = 10.0f - nextafterf(10, 0);
    for (float x = 10.0f; x >= 0.0f; x -= dx)
    {
        N++;

        Vec3 r_ij = {
            .x = x,
            .y = x,
            .z = x,
        };

        // First image convention
        r_ij.x = minimum_image(r_ij.x, 1);
        r_ij.y = minimum_image(r_ij.y, 1);
        r_ij.z = minimum_image(r_ij.z, 1);

        double r_ij_mod = r_ij.x * r_ij.x + r_ij.y * r_ij.y + r_ij.z * r_ij.z;
    }
    c1 = clock();
    printf("%d time: %g\n", N, (double)(c1 - c0) / CLOCKS_PER_SEC);
}