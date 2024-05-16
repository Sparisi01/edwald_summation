#include "constants.h"
#include "structures.h"
#include <math.h>
#include <stdlib.h>

#ifndef THERMODINAMICS_H
#define THERMODINAMICS_H

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

// The first time radialDistribution is called
// alloc and array of bins, each bin correspond to a certain radial range.
// Each time the function is called calculate the radial distribution
// of particle around each other. For each couple of particle the distance is calculated
// and the corresponding bin is increased by 1. The bin array is independent by the single
// function call. This let you make statistic on different system time.
// To reset the array pass the OPTION = 0.
double *radialDistribution(System *system, short OPTION)
{
    static int *bins;
    static int has_been_allocated = 0;

    const int n_bins = 1;
    const double rad_step = CELL_LENGHT / n_bins;

    if (!has_been_allocated)
    {
        bins = (int *)calloc(n_bins, sizeof(int));
        if (!bins)
        {
            perror("Error allocating bins array");
            exit(EXIT_FAILURE);
        }
    }

    switch (OPTION)
    {
    case 0: // reset the array and then do case 1
        free(bins);
        bins = (int *)calloc(n_bins, sizeof(int));
        if (!bins)
        {
            perror("Error allocating bins array");
            exit(EXIT_FAILURE);
        }
    case 1: // add current system radial distribution to the bins array
        // TODO: perform radial distribution
        break;
    default:
        break;
    }

    return 0;
}

#endif