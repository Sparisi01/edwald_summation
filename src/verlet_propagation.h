#include "structures.h"
#include <stdlib.h>

#ifndef VERLET_PROPAGATION_H
#define VERLET_PROPAGATION_H

int verletPropagationStep(System *system, double time_step, Vec3 *(*forceFunction)(System *system, double *args), double *args)
{

    for (size_t i = 0; i < system->n_particles; i++)
    {
        system->particles[i].x += system->particles[i].vx * time_step + 0.5 / system->particles[i].mass * system->forces[i].x * time_step * time_step;
        system->particles[i].y += system->particles[i].vy * time_step + 0.5 / system->particles[i].mass * system->forces[i].y * time_step * time_step;
        system->particles[i].z += system->particles[i].vz * time_step + 0.5 / system->particles[i].mass * system->forces[i].z * time_step * time_step;
    }

    system->time += time_step;

    Vec3 *new_forces = forceFunction(system, args);
    if (!new_forces)
    {
        return 1;
    }

    for (size_t i = 0; i < system->n_particles; i++)
    {
        system->particles[i].vx += 0.5 / system->particles[i].mass * (system->forces[i].x + new_forces[i].x) * time_step;
        system->particles[i].vy += 0.5 / system->particles[i].mass * (system->forces[i].y + new_forces[i].y) * time_step;
        system->particles[i].vz += 0.5 / system->particles[i].mass * (system->forces[i].z + new_forces[i].z) * time_step;
    }

    free(system->forces);
    system->forces = new_forces;

    return 0;
}

#endif