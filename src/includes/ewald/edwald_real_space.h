#if !defined(EDWALD_REal_SPACE_H)
#define EDWALD_REal_SPACE_H

#include "../constants.h"
#include "../structures.h"
#include <math.h>
#include <stdlib.h>

double _CUTOFF = NAN;

double real_space_potential_coulomb(System *s, double ALPHA)
{
    double sum = 0;
    double qi, qj;

    // Using the fact that Uᵢⱼ = Uⱼᵢ sum over i and j with j>i and remove the 1/2 factor.
    // By doing that you cut in half the number of CPU cicles needed
    for (size_t i = 0; i < s->n_particles - 1; i++)
    {
        for (size_t j = i + 1; j < s->n_particles; j++)
        {
            qi = s->particles[i].charge;
            qj = s->particles[j].charge;

            Vec3 r_ij = {
                .x = s->particles[i].x - s->particles[j].x,
                .y = s->particles[i].y - s->particles[j].y,
                .z = s->particles[i].z - s->particles[j].z,
            };

            // First image convention
            r_ij.x = r_ij.x - rint((r_ij.x) / s->cell_lenght) * s->cell_lenght;
            r_ij.y = r_ij.y - rint((r_ij.y) / s->cell_lenght) * s->cell_lenght;
            r_ij.z = r_ij.z - rint((r_ij.z) / s->cell_lenght) * s->cell_lenght;

            double r_ij_mod = sqrt(r_ij.x * r_ij.x + r_ij.y * r_ij.y + r_ij.z * r_ij.z);

            // In first image convention _CUTOFF must be less than L/2
            if (r_ij_mod > _CUTOFF) continue;

            sum += (qi * qj) * erfc(ALPHA * r_ij_mod) / r_ij_mod;
        }
    }
    return sum;
}

double real_space_potential_yukawa(System *s, double ALPHA, double LAMBDA)
{
    double sum = 0;
    double qi, qj;

    // Using the fact that Uᵢⱼ = Uⱼᵢ sum over i and j with j>i and remove the 1/2 factor.
    // By doing that you cut in half the number of CPU cicles needed
    for (size_t i = 0; i < s->n_particles - 1; i++)
    {
        for (size_t j = i + 1; j < s->n_particles; j++)
        {
            qi = s->particles[i].charge;
            qj = s->particles[j].charge;

            Vec3 r_ij = {
                .x = s->particles[i].x - s->particles[j].x,
                .y = s->particles[i].y - s->particles[j].y,
                .z = s->particles[i].z - s->particles[j].z,
            };

            // First image convention
            r_ij.x = r_ij.x - rint((r_ij.x) / s->cell_lenght) * s->cell_lenght;
            r_ij.y = r_ij.y - rint((r_ij.y) / s->cell_lenght) * s->cell_lenght;
            r_ij.z = r_ij.z - rint((r_ij.z) / s->cell_lenght) * s->cell_lenght;

            double r_ij_mod = sqrt(r_ij.x * r_ij.x + r_ij.y * r_ij.y + r_ij.z * r_ij.z);

            // In first image convention _CUTOFF must be less than L/2
            if (r_ij_mod > _CUTOFF) continue;

            double t1 = erfc(ALPHA * r_ij_mod + 0.5 * LAMBDA / ALPHA) * exp(LAMBDA * r_ij_mod);
            double t2 = erfc(ALPHA * r_ij_mod - 0.5 * LAMBDA / ALPHA) * exp(-LAMBDA * r_ij_mod);
            sum += (qi * qj) * (t1 + t2) / (2 * r_ij_mod);
        }
    }
    return sum;
}

Vec3 *real_space_force(System *s)
{
}

#endif // EDWALD_REal_SPACE_H
