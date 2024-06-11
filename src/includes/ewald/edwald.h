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

#endif