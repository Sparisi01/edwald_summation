#if !defined(EDWALD_REal_SPACE_H)
#define EDWALD_REal_SPACE_H

#include "../constants.h"
#include "../structures.h"
#include <math.h>
#include <stdlib.h>

double _CUTOFF = INFINITY;

double real_space_potential_coulomb(System *s, double ALPHA)
{
    double sum = 0;
    for (size_t i = 0; i < s->n_particles - 1; i++)
    {
        for (size_t j = i + 1; j < s->n_particles; j++)
        {
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

            if (r_ij_mod > _CUTOFF) continue;

            sum += (s->particles[i].charge * s->particles[j].charge) * erfc(ALPHA * r_ij_mod) / r_ij_mod;
        }
    }
    return sum;
}

Vec3 *real_space_force(System *s)
{
}

#endif // EDWALD_REal_SPACE_H
