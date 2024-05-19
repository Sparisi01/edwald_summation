#if !defined(COULMB_POTENTIAL_H)
#define COULMB_POTENTIAL_H

#include <math.h>
#include <stdlib.h>

#include "../constants.h"
#include "../structures.h"

int _R_RANGE = 5;

double getCoulombPotential(System *s)
{
    double sum = 0;

    for (size_t i = 0; i < s->n_particles; i++)
    {
        for (size_t j = i; j < s->n_particles; j++)
        {
            for (int r_x = -_R_RANGE; r_x <= _R_RANGE; r_x++)
            {
                for (int r_y = -_R_RANGE; r_y <= _R_RANGE; r_y++)
                {
                    for (int r_z = -_R_RANGE; r_z <= _R_RANGE; r_z++)
                    {
                        // Exclude self particle in cell (0,0,0)
                        if (r_x == 0 && r_y == 0 && r_z == 0 && i == j) continue;

                        Vec3 r_ij = {
                            .x = s->particles[i].x - s->particles[j].x,
                            .y = s->particles[i].y - s->particles[j].y,
                            .z = s->particles[i].z - s->particles[j].z,
                        };

                        r_ij.x = r_ij.x - rint((r_ij.x) / s->cell_lenght) * s->cell_lenght + r_x * s->cell_lenght;
                        r_ij.y = r_ij.y - rint((r_ij.y) / s->cell_lenght) * s->cell_lenght + r_y * s->cell_lenght;
                        r_ij.z = r_ij.z - rint((r_ij.z) / s->cell_lenght) * s->cell_lenght + r_z * s->cell_lenght;

                        double r_ij_mod = sqrt(r_ij.x * r_ij.x + r_ij.y * r_ij.y + r_ij.z * r_ij.z);
                        sum += (s->particles[i].charge * s->particles[j].charge) / r_ij_mod;
                    }
                }
            }
        }
    }

    return sum * POT_TYPE;
}

double getYukawaPotential(System *s, double lambda)
{
    double sum = 0;

    for (size_t i = 0; i < s->n_particles; i++)
    {
        for (size_t j = i; j < s->n_particles; j++)
        {
            for (int r_x = -_R_RANGE; r_x <= _R_RANGE; r_x++)
            {
                for (int r_y = -_R_RANGE; r_y <= _R_RANGE; r_y++)
                {
                    for (int r_z = -_R_RANGE; r_z <= _R_RANGE; r_z++)
                    {
                        // Exclude self particle in cell (0,0,0)
                        if (r_x == 0 && r_y == 0 && r_z == 0 && i == j) continue;

                        Vec3 r_ij = {
                            .x = s->particles[i].x - s->particles[j].x,
                            .y = s->particles[i].y - s->particles[j].y,
                            .z = s->particles[i].z - s->particles[j].z,
                        };

                        r_ij.x = r_ij.x - rint((r_ij.x) / s->cell_lenght) * s->cell_lenght + r_x * s->cell_lenght;
                        r_ij.y = r_ij.y - rint((r_ij.y) / s->cell_lenght) * s->cell_lenght + r_y * s->cell_lenght;
                        r_ij.z = r_ij.z - rint((r_ij.z) / s->cell_lenght) * s->cell_lenght + r_z * s->cell_lenght;

                        double r_ij_mod = sqrt(r_ij.x * r_ij.x + r_ij.y * r_ij.y + r_ij.z * r_ij.z);

                        sum += (s->particles[i].charge * s->particles[j].charge) / r_ij_mod * exp(-lambda * r_ij_mod);
                    }
                }
            }
        }
    }

    return sum * POT_TYPE;
}

#endif // COULMB_POTENTIAL_H
