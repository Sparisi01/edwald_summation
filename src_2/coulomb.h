#ifndef COULOMB_POTENTIAL_H
#define COULOMB_POTENTIAL_H

#include "constants.h"
#include "structures.h"

#include <math.h>
#include <stdlib.h>

int COULOMB_CELL_RANGE = 3;

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

    for (size_t i = 0; i < n_particles; i++)
    {
        for (size_t j = i; j < n_particles; j++)
        {
            for (int n_k_x = -COULOMB_CELL_RANGE; n_k_x <= COULOMB_CELL_RANGE; n_k_x++)
            {
                for (int n_k_y = -COULOMB_CELL_RANGE; n_k_y <= COULOMB_CELL_RANGE; n_k_y++)
                {
                    for (int n_k_z = -COULOMB_CELL_RANGE; n_k_z <= COULOMB_CELL_RANGE; n_k_z++)
                    {

                        // Particle i == j must be removed only in the (0,0,0) cell
                        if (i == j && n_k_x == 0 && n_k_y == 0 && n_k_z == 0) continue;

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

                        // NO third principle if the particle is the same
                        if (i == j) continue;

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

double coulomb_potential(System *system)
{

    int n_particles = system->n_particles;
    Particle *particles = system->particles;

    double pot_energy = 0;
    for (size_t i = 0; i < n_particles - 1; i++)
    {
        for (size_t j = i; j < n_particles; j++)
        {
            for (int n_k_x = -COULOMB_CELL_RANGE; n_k_x <= COULOMB_CELL_RANGE; n_k_x++)
            {
                for (int n_k_y = -COULOMB_CELL_RANGE; n_k_y <= COULOMB_CELL_RANGE; n_k_y++)
                {
                    for (int n_k_z = -COULOMB_CELL_RANGE; n_k_z <= COULOMB_CELL_RANGE; n_k_z++)
                    {
                        if (i == j && n_k_x == 0 && n_k_y == 0 && n_k_z == 0) continue;

                        double r_ij_x = particles[i].x - particles[j].x - rint((particles[i].x - particles[j].x) / CELL_LENGHT) * CELL_LENGHT + n_k_x * CELL_LENGHT; //
                        double r_ij_y = particles[i].y - particles[j].y - rint((particles[i].y - particles[j].y) / CELL_LENGHT) * CELL_LENGHT + n_k_y * CELL_LENGHT; //
                        double r_ij_z = particles[i].z - particles[j].z - rint((particles[i].z - particles[j].z) / CELL_LENGHT) * CELL_LENGHT + n_k_z * CELL_LENGHT; //
                        double r_ij = sqrt(r_ij_x * r_ij_x + r_ij_y * r_ij_y + r_ij_z * r_ij_z);

                        pot_energy += FORCE_TYPE_CONSTANT * particles[i].charge * particles[j].charge / r_ij * exp(-LAMBDA * r_ij);
                    }
                }
            }
        }
    }
    return pot_energy;
}

#endif