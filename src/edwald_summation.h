#include "constants.h"
#include "structures.h"
#include <math.h>
#include <stdlib.h>

#ifndef EDWALD_SUMMATION_H
#define EDWALD_SUMMATION_H

int restore_positions_in_lattice_first_cell(Particle *particles, int n_particles)
{
    for (size_t i = 0; i < n_particles; i++)
    {
        particles[i].x = particles[i].x - rint(particles[i].x / CELL_LENGHT) * CELL_LENGHT;
        particles[i].y = particles[i].y - rint(particles[i].y / CELL_LENGHT) * CELL_LENGHT;
        particles[i].z = particles[i].z - rint(particles[i].z / CELL_LENGHT) * CELL_LENGHT;
    }
    return 0;
}

Vec3 *edwald_summation(System *system, double *args)
{
    Particle *particles = system->particles;

    Vec3 *tmp_forces = (Vec3 *)malloc(system->n_particles * sizeof(Vec3));
    if (!tmp_forces)
    {
        return NULL;
    }

    for (size_t i = 0; i < system->n_particles; i++)
    {
        tmp_forces[i].x = 0;
        tmp_forces[i].y = 0;
        tmp_forces[i].z = 0;
    }

    // Bring back all the particles in the (0,0,0) cell
    restore_positions_in_lattice_first_cell(particles, system->n_particles);

    const int SELECTED_PARTICLE = 0;
    static int r_k_print_countdown = 5;
    double real_sum = 0;
    double k_sum = 0;

    for (size_t i = 0; i < system->n_particles - 1; i++)
    {
        for (size_t j = i + 1; j < system->n_particles; j++)
        {
            double r_ij_x = particles[i].x - particles[j].x - rint((particles[i].x - particles[j].x) / CELL_LENGHT) * CELL_LENGHT; //
            double r_ij_y = particles[i].y - particles[j].y - rint((particles[i].y - particles[j].y) / CELL_LENGHT) * CELL_LENGHT; //
            double r_ij_z = particles[i].z - particles[j].z - rint((particles[i].z - particles[j].z) / CELL_LENGHT) * CELL_LENGHT; //
            double r_ij = sqrt(r_ij_x * r_ij_x + r_ij_y * r_ij_y + r_ij_z * r_ij_z);                                               // |r_i - r_j|

            /* if (r_ij >= CELL_L / 2) {
                continue;  // CUTOF potential
            } */

            double classical_force = particles[i].charge * particles[j].charge / (r_ij * r_ij * r_ij);                         // Classic Coulomb Like force
            double edwald_correction = 2 * ALPHA / SQR_PI * exp(-(ALPHA * r_ij) * (ALPHA * r_ij)) * r_ij + erfc(ALPHA * r_ij); // Edwald correction in real space

            // Force on particle i due to j
            tmp_forces[i].x += FORCE_TYPE_CONSTANT * r_ij_x * classical_force * edwald_correction;
            tmp_forces[i].y += FORCE_TYPE_CONSTANT * r_ij_y * classical_force * edwald_correction;
            tmp_forces[i].z += FORCE_TYPE_CONSTANT * r_ij_z * classical_force * edwald_correction;

            // tmp_forces on particle j due to i using third principle
            tmp_forces[j].x -= FORCE_TYPE_CONSTANT * r_ij_x * classical_force * edwald_correction;
            tmp_forces[j].y -= FORCE_TYPE_CONSTANT * r_ij_y * classical_force * edwald_correction;
            tmp_forces[j].z -= FORCE_TYPE_CONSTANT * r_ij_z * classical_force * edwald_correction;

            if (i == SELECTED_PARTICLE)
            {
                real_sum += FORCE_TYPE_CONSTANT * r_ij_x * classical_force * edwald_correction;
                real_sum += FORCE_TYPE_CONSTANT * r_ij_y * classical_force * edwald_correction;
                real_sum += FORCE_TYPE_CONSTANT * r_ij_z * classical_force * edwald_correction;
            }
        }

        // SECTION - Reciprocal space force

        /**NOTE - La formula non prevede l'esclusione dell'interazione tra particella i e j
         * se esse appartengono a celle diverse. Noto però che in questo caso il sin da valore
         * nullo. Posso escludere l'interazione.
         **/
        for (size_t j = i + 1; j < system->n_particles; j++)
        {
            for (int n_k_x = -N_K_RANGE; n_k_x <= N_K_RANGE; n_k_x++)
            {
                for (int n_k_y = -N_K_RANGE; n_k_y <= N_K_RANGE; n_k_y++)
                {
                    for (int n_k_z = -N_K_RANGE; n_k_z <= N_K_RANGE; n_k_z++)
                    {
                        if (n_k_x == 0 && n_k_y == 0 && n_k_z == 0)
                        {
                            continue;
                        }
                        double k_x = n_k_x * 2 * PI / CELL_LENGHT;               // Componente x vettore spazione reciproco
                        double k_y = n_k_y * 2 * PI / CELL_LENGHT;               // Componente y vettore spazione reciproco
                        double k_z = n_k_z * 2 * PI / CELL_LENGHT;               // Componente z vettore spazione reciproco
                        double k_vec_mag2 = (k_x * k_x + k_y * k_y + k_z * k_z); // Magnitude square of reciprocal-lattice vector

                        double dot_prod_pos_k = (k_x * (system->particles[i].x - system->particles[j].x) + k_y * (system->particles[i].y - system->particles[j].y) + k_z * (system->particles[i].z - system->particles[j].z));
                        double long_distance_force = 4 * PI / CELL_VOLUME * system->particles[i].charge * system->particles[j].charge;
                        long_distance_force *= exp(-k_vec_mag2 / (2 * ALPHA * ALPHA)) / k_vec_mag2 * sin(dot_prod_pos_k);

                        // tmp_forces on particle i due to j copies
                        tmp_forces[i].x += FORCE_TYPE_CONSTANT * k_x * long_distance_force;
                        tmp_forces[i].z += FORCE_TYPE_CONSTANT * k_y * long_distance_force;
                        tmp_forces[i].y += FORCE_TYPE_CONSTANT * k_z * long_distance_force;
                        // tmp_forces on particle j due to i copies using third principle
                        tmp_forces[j].x -= FORCE_TYPE_CONSTANT * k_x * long_distance_force;
                        tmp_forces[j].y -= FORCE_TYPE_CONSTANT * k_y * long_distance_force;
                        tmp_forces[j].z -= FORCE_TYPE_CONSTANT * k_z * long_distance_force;

                        if (i == SELECTED_PARTICLE)
                        {
                            k_sum += FORCE_TYPE_CONSTANT * k_x * long_distance_force;
                            k_sum += FORCE_TYPE_CONSTANT * k_y * long_distance_force;
                            k_sum += FORCE_TYPE_CONSTANT * k_z * long_distance_force;
                        }
                    }
                }
            }
        }
        /*!SECTION */
    }
    if (r_k_print_countdown > 0)
    {
        printf("REAL SUM: %.10E\nRECIPROCAL SUM: %.10E\n", real_sum, k_sum);
        r_k_print_countdown--;
    }

    return tmp_forces;
}
#endif