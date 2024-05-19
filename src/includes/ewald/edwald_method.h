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

double self_potential(System *system, double ALPHA)
{
    double sum = 0;
    for (size_t i = 0; i < system->n_particles; i++)
    {
        sum += pow(system->particles[i].charge, 2);
    }
    return sum * ALPHA / SQR_PI;
}

double getEdwaldPotential(System *system)
{
    double short_range = real_space_potential(system, _ALPHA);
    double long_range = reciprocal_space_potential(system, _ALPHA);
    double self = self_potential(system, _ALPHA);

    return short_range + long_range - self;
}

Vec3 *getEdwaldForces(System *system)
{
}

#endif