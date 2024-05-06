#include "structures.h"

#ifndef OBSERVABLES_H
#define OBSERVABLES_H

double kinetic_energy(Particle *particles, int n_particles)
{
    double kinetic_energy = 0;

    for (size_t i = 0; i < n_particles; i++)
    {
        kinetic_energy += 0.5 * particles[i].mass * (particles[i].vx * particles[i].vx);
        kinetic_energy += 0.5 * particles[i].mass * (particles[i].vy * particles[i].vy);
        kinetic_energy += 0.5 * particles[i].mass * (particles[i].vz * particles[i].vz);
    }

    return kinetic_energy;
}

double potential_energy(Particle *particles, int n_particles)
{
    double pot_energy = 0;
    for (size_t i = 0; i < n_particles - 1; i++)
    {
        for (size_t j = i + 1; j < n_particles; j++)
        {
            double r_ij_x = particles[i].x - particles[j].x - rint((particles[i].x - particles[j].x) / CELL_LENGHT) * CELL_LENGHT; //
            double r_ij_y = particles[i].y - particles[j].y - rint((particles[i].y - particles[j].y) / CELL_LENGHT) * CELL_LENGHT; //
            double r_ij_z = particles[i].z - particles[j].z - rint((particles[i].z - particles[j].z) / CELL_LENGHT) * CELL_LENGHT; //
            double r_ij = sqrt(r_ij_x * r_ij_x + r_ij_y * r_ij_y + r_ij_z * r_ij_z);

            pot_energy += FORCE_TYPE_CONSTANT * particles[i].charge * particles[j].charge / r_ij;
        }
    }
    return pot_energy;
}

#endif