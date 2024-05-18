#include "edwald_real_space.h"
#include "edwald_reciprocal_space.h"
#include <math.h>

#ifndef EDWALD_H
#define EDWALD_H

double _ALPHA = 1;

void setAlphaByError(double error)
{
    _ALPHA = sqrt(-log(error));
}

double getEdwaldPotential(System *system)
{
}

Vec3 *getEdwaldForces(System *system)
{
}

#endif