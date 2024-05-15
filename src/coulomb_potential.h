#ifndef COULOMB_POTENTIAL_H
#define COULOMB_POTENTIAL_H

#include "constants.h"
#include "structures.h"

#include <math.h>
#include <stdlib.h>

#define CELLS_RANGE 0

Vec3 *coulomb_force(System *system, double *args)
{

    Particle *particles = system->particles;
    int n_particles = system->n_particles;

    Vec3 *tmp_forces = (Vec3 *)malloc(n_particles * sizeof(Vec3));
    if (!tmp_forces)
    {
        return NULL;
    }

    for (size_t i = 0; i < n_particles; i++)
    {
        tmp_forces[i].x = 0;
        tmp_forces[i].y = 0;
        tmp_forces[i].z = 0;
    }

    for (size_t i = 0; i < n_particles - 1; i++)
    {
        for (size_t j = i + 1; j < n_particles; j++)
        {

            for (int n_k_x = -CELLS_RANGE; n_k_x <= CELLS_RANGE; n_k_x++)
            {
                for (int n_k_y = -CELLS_RANGE; n_k_y <= CELLS_RANGE; n_k_y++)
                {
                    for (int n_k_z = -CELLS_RANGE; n_k_z <= CELLS_RANGE; n_k_z++)
                    {
                        Vec3 r_ij = {
                            .x = particles[i].x - particles[j].x,
                            .y = particles[i].y - particles[j].y,
                            .z = particles[i].z - particles[j].z};

                        r_ij.x = r_ij.x - rint((r_ij.x) / CELL_LENGHT) * CELL_LENGHT + n_k_x * CELL_LENGHT;
                        r_ij.y = r_ij.y - rint((r_ij.y) / CELL_LENGHT) * CELL_LENGHT + n_k_y * CELL_LENGHT;
                        r_ij.z = r_ij.z - rint((r_ij.z) / CELL_LENGHT) * CELL_LENGHT + n_k_z * CELL_LENGHT;

                        double r_ij_mod = sqrt(r_ij.x * r_ij.x + r_ij.y * r_ij.y + r_ij.z * r_ij.z);
                        double mod3 = r_ij_mod * r_ij_mod * r_ij_mod;
                        Vec3 force_ij = {.x = r_ij.x / mod3, .y = r_ij.y / mod3, .z = r_ij.z / mod3};

                        tmp_forces[i].x += FORCE_TYPE_CONSTANT * particles[i].charge * particles[j].charge * force_ij.x;
                        tmp_forces[i].y += FORCE_TYPE_CONSTANT * particles[i].charge * particles[j].charge * force_ij.y;
                        tmp_forces[i].z += FORCE_TYPE_CONSTANT * particles[i].charge * particles[j].charge * force_ij.z;

                        // Use Newton third principle
                        tmp_forces[j].x -= FORCE_TYPE_CONSTANT * particles[i].charge * particles[j].charge * force_ij.x;
                        tmp_forces[j].y -= FORCE_TYPE_CONSTANT * particles[i].charge * particles[j].charge * force_ij.y;
                        tmp_forces[j].z -= FORCE_TYPE_CONSTANT * particles[i].charge * particles[j].charge * force_ij.z;
                    }
                }
            }
        }
    }
    return tmp_forces;
}

#endif