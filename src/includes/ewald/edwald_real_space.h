#if !defined(EDWALD_REal_SPACE_H)
#define EDWALD_REal_SPACE_H

#include "../constants.h"
#include "../structures.h"
#include "../utils/lattice.h"
#include <inttypes.h>
#include <math.h>
#include <stdlib.h>

double _CUTOFF = 10;
int _R_RANGE_EWALD = 0;
int _VERBOSE = 0;
// int _R_RANGE_EWALD = 0; // When set to 0 first image convention

// _R_RANGE_EWALD defines the rnage in real space, if _CUTOFF not  set to infinity it is applied
// and the sum is on spheres and not squares.
// There is no removal of double interaction thus the final energy must be divided
// by 2.

double real_space_coulomb_energy(System *s, double ALPHA)
{
    _R_RANGE_EWALD = floor(_CUTOFF / (s->cell_lenght / 2)); // Set real range based on cutoff in real space
    // int _R_RANGE_EWALD = 10;
    if (ALPHA <= 0) goto ALPHA_ERROR;
    if (_CUTOFF <= 0) goto CUTOFF_ERROR;

    int64_t skipped_for_truncation = 0;
    double energy_sum = 0;
    double qi, qj;

    // Only cell (0,0,0), cutoff must be less than cell_length
    for (size_t i = 0; i < s->n_particles; i++)
    {
        qi = s->particles[i].charge;

        for (size_t j = 0; j < s->n_particles; j++)
        {
            for (int r_x = -_R_RANGE_EWALD; r_x <= _R_RANGE_EWALD; r_x++)
            {
                for (int r_y = -_R_RANGE_EWALD; r_y <= _R_RANGE_EWALD; r_y++)
                {
                    for (int r_z = -_R_RANGE_EWALD; r_z <= _R_RANGE_EWALD; r_z++)
                    {
                        // Exclude self particle in cell (0,0,0)
                        if (r_x == 0 && r_y == 0 && r_z == 0 && i == j) continue;
                        // if (i == j) continue;

                        qj = s->particles[j].charge;

                        Vec3 r_ij = {
                            .x = s->particles[i].x - s->particles[j].x,
                            .y = s->particles[i].y - s->particles[j].y,
                            .z = s->particles[i].z - s->particles[j].z,
                        };

                        // First image convention
                        r_ij.x = minimum_image(r_ij.x, s->cell_lenght) + r_x * s->cell_lenght;
                        r_ij.y = minimum_image(r_ij.y, s->cell_lenght) + r_y * s->cell_lenght;
                        r_ij.z = minimum_image(r_ij.z, s->cell_lenght) + r_z * s->cell_lenght;

                        double r_ij_mod = sqrt(r_ij.x * r_ij.x + r_ij.y * r_ij.y + r_ij.z * r_ij.z);

                        // In first image convention _CUTOFF must be less than L/2
                        if (r_ij_mod > _CUTOFF) continue;
                        if (r_ij_mod == 0) goto PARTICLE_OVERLAP_ERROR;

                        double old_energy_sum = energy_sum;
                        energy_sum += (qi * qj) * erfc(ALPHA * r_ij_mod) / r_ij_mod;

                        if (old_energy_sum == energy_sum && qi * qj != 0)
                        {
                            skipped_for_truncation += 1;
                        }
                    }
                }
            }
        }
    }

    if (skipped_for_truncation != 0 && _VERBOSE == 1) printf("WARNING: %lld skipped for truncation\n", skipped_for_truncation);
    return 0.5 * energy_sum;

PARTICLE_OVERLAP_ERROR:
    printf("ERROR: particle overlap");
    exit(EXIT_FAILURE);

ALPHA_ERROR:
    printf("ERROR: ALPHA must be greater than 0");
    exit(EXIT_FAILURE);

CUTOFF_ERROR:
    printf("ERROR: _CUTOFF must be greater than 0");
    exit(EXIT_FAILURE);
}

Vec3 *real_space_force(System *s)
{
    return NULL;
}

#endif // EDWALD_REal_SPACE_H
