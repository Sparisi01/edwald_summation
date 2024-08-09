#ifndef EDWALD_H
#define EDWALD_H

#include "edwald_real_space.h"
#include "edwald_reciprocal_space.h"
#include <math.h>

double _ALPHA = 1;

double alpha_by_precision(double precision)
{
    return sqrt(-log(precision));
}

double self_coulomb_energy(System *system, double ALPHA)
{
    double sum = 0;
    for (size_t i = 0; i < system->n_particles; i++)
    {
        sum += pow(system->particles[i].charge, 2);
    }
    return sum * ALPHA / SQR_PI;
}

double ewald_energy(System *system)
{
    double short_range = real_space_coulomb_energy(system, _ALPHA);
    double long_range = reciprocal_space_coulomb_energy(system, _ALPHA);
    double self = self_coulomb_energy(system, _ALPHA);

    return (short_range + long_range - self);
}

Vec3 *getEdwaldForces(System *system)
{
    return NULL;
}

double errorsDifference(double error, double s, double Q, double cell_length, double alpha)
{
    return exp(-(s * s)) / (pow(s, 3. / 2)) * Q * sqrt((2 + PI) / (2 * PI * alpha * pow(cell_length, 3))) - error;
};

double findSbybisection(double a, double b, double error, double Q, double cell_length, double alpha, double precision)
{
    double max = errorsDifference(error, a, Q, cell_length, alpha);
    double min = errorsDifference(error, b, Q, cell_length, alpha);

    if (!(max > 0 && min < 0))
    {
        printf("Max e min non rispettano parametri bisezione\n");
    };
    double c = 0;
    int root_find = 0;
    while (!root_find)
    {
        c = (a + b) / 2;
        double mid_value = errorsDifference(error, c, Q, cell_length, alpha);
        // printf("%.5E\n", mid_value);
        if (fabs(mid_value) < precision)
        {
            root_find = 1;
        }
        else
        {
            if (mid_value < 0) b = c;
            if (mid_value > 0) a = c;
        }
    }

    printf("ROOT FOUND: %.5E\n", c);
    return c;
}

double optimizeParameter(double error, double cell_length, int N_particles, double Q)
{
    double alpha = pow((TAU_RAPP * pow(PI, 3) / pow(cell_length, 6) * N_particles), 1. / 6);

    printf("ALPHA: %.5E\n", alpha);

    double s = findSbybisection(1e-3, 1e3, error, Q, cell_length, alpha, 1e-10);

    double rc = s / alpha;
    double kc = 2 * s * alpha;

    _ALPHA = alpha;
    _CUTOFF = rc;
    _K_RANGE_EWALD = ceil(kc / (2 * PI / cell_length));

    printf("OPTIMIZED PARAMETERS: R_C = %.5E, N_C = %.5E\n", rc, kc / (2 * PI / cell_length));
}

#endif