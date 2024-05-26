#include "edwald_real_space.h"
#include "edwald_reciprocal_space.h"
#include <math.h>

#ifndef EDWALD_H
#define EDWALD_H

double _ALPHA = 1; // ciao

void setAlphaByError(double error)
{
    _ALPHA = sqrt(-log(error));
}

double self_potential_coulomb(System *system, double ALPHA)
{
    double sum = 0;
    for (size_t i = 0; i < system->n_particles; i++)
    {
        sum += pow(system->particles[i].charge, 2);
    }
    return sum * ALPHA / SQR_PI;
}

double getEdwaldPotentialCoulomb(System *system)
{
    double short_range = real_space_potential_coulomb(system, _ALPHA);
    double long_range = reciprocal_space_potential_coulomb(system, _ALPHA);
    double self = self_potential_coulomb(system, _ALPHA);

    return (short_range + long_range - self) * POT_TYPE;
}

double getEdwaldPotentialYukawa(System *system, double LAMBDA)
{
    double short_range = real_space_potential_yukawa(system, _ALPHA, LAMBDA);
    double long_range = reciprocal_space_potential_yukawa(system, _ALPHA, LAMBDA);

    return (short_range + long_range) * POT_TYPE;
}

Vec3 *getEdwaldForces(System *system)
{
}

#endif