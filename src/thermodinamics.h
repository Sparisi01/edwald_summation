#include <math.h>

#include "constants.h"

double temperature(double *vel)
{
    return 1;
}

double kinetic_energy(double *vel, double *masses, int n_particles)
{
    double kinetic_energy = 0;

    for (size_t i = 0; i < n_particles * 3; i += 3)
    {
        kinetic_energy += 0.5 * masses[i / 3] * (vel[i + 0] * vel[i + 0] + vel[i + 1] * vel[i + 1] + vel[i + 2] * vel[i + 2]);
    }

    return kinetic_energy;
}

double potential_energy(double *pos, double *particle_charge, int n_particles)
{
    double pot_energy = 0;
    for (size_t i = 0; i < n_particles * 3; i += 3)
    {
        for (size_t j = i + 3; j < n_particles * 3; j += 3)
        {
            double r_ij_x = pos[i + 0] - pos[j + 0] - rint((pos[i + 0] - pos[j + 0]) / CELL_L) * CELL_L; //
            double r_ij_y = pos[i + 1] - pos[j + 1] - rint((pos[i + 1] - pos[j + 1]) / CELL_L) * CELL_L; //
            double r_ij_z = pos[i + 2] - pos[j + 2] - rint((pos[i + 2] - pos[j + 2]) / CELL_L) * CELL_L; //
            double vec_dist_mag = sqrt(pow(r_ij_x, 2) + pow(r_ij_y, 2) + pow(r_ij_z, 2));
            pot_energy += FORCE_TYPE_CONSTANT * particle_charge[i / 3] * particle_charge[j / 3] / vec_dist_mag;
        }
    }
    return pot_energy;
}
