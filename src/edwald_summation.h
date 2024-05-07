#include "constants.h"
#include "interpolation.h"
#include "structures.h"
#include <math.h>
#include <stdlib.h>
#include <time.h>

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

Vec3 compute_reciprocal_space_force(double r_ij_x, double r_ij_y, double r_ij_z)
{
    Vec3 force = {.x = 0, .y = 0, .z = 0};

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

                double dot_prod_pos_k = (k_x * r_ij_x + k_y * r_ij_y + k_z * r_ij_z);
                double long_distance_force = 4 * PI / CELL_VOLUME;
                long_distance_force *= exp(-k_vec_mag2 / (2 * ALPHA * ALPHA)) / k_vec_mag2 * sin(dot_prod_pos_k);

                // tmp_forces on particle i due to j copies
                force.x += FORCE_TYPE_CONSTANT * k_x * long_distance_force;
                force.z += FORCE_TYPE_CONSTANT * k_y * long_distance_force;
                force.y += FORCE_TYPE_CONSTANT * k_z * long_distance_force;
            }
        }
    }

    return force;
}

Vec3 compute_real_space_force(double r_ij_x, double r_ij_y, double r_ij_z)
{
    Vec3 force = {.x = 0, .y = 0, .z = 0};

    // In real space i-force has to be compute with the j-particle
    // nearest copy in the lattice.
    r_ij_x = r_ij_x - rint((r_ij_x) / CELL_LENGHT) * CELL_LENGHT;
    r_ij_y = r_ij_y - rint((r_ij_y) / CELL_LENGHT) * CELL_LENGHT;
    r_ij_z = r_ij_z - rint((r_ij_z) / CELL_LENGHT) * CELL_LENGHT;

    double r_ij = sqrt(r_ij_x * r_ij_x + r_ij_y * r_ij_y + r_ij_z * r_ij_z);

    /* if (r_ij >= CELL_L / 2) {
        continue;  // CUTOF potential
    } */

    double edwald_correction = 2 * ALPHA / SQR_PI * exp(-(ALPHA * r_ij) * (ALPHA * r_ij)) * r_ij + erfc(ALPHA * r_ij); // Edwald correction in real space

    force.x = r_ij_x / (r_ij * r_ij * r_ij) * edwald_correction;
    force.y = r_ij_y / (r_ij * r_ij * r_ij) * edwald_correction;
    force.z = r_ij_z / (r_ij * r_ij * r_ij) * edwald_correction;

    return force;
}

Vec3 tabulated_reciprocal_space_term(double r_ij_x, double r_ij_y, double r_ij_z)
{
    const int table_size = RECIPROCAL_SPACE_TABLE_SIZE;
    const double table_step = 2 * CELL_LENGHT / RECIPROCAL_SPACE_TABLE_SIZE;
    static int has_been_tabled = 0;
    static Vec3 table
        [RECIPROCAL_SPACE_TABLE_SIZE + 4]
        [RECIPROCAL_SPACE_TABLE_SIZE + 4]
        [RECIPROCAL_SPACE_TABLE_SIZE + 4];

    if (!has_been_tabled)
    {
        printf("STARTING GENERATION TABLE\n");
        for (int i = 0; i < table_size + 4; i++)
        {
            for (int j = 0; j < table_size; j++)
            {
                for (int k = 0; k < table_size; k++)
                {
                    table[i][j][k] = compute_reciprocal_space_force(
                        table_step * (i - table_size / 2),
                        table_step * (j - table_size / 2),
                        table_step * (k - table_size / 2));
                }
            }
        }
        has_been_tabled = 1;
        printf("TABLE COMPLETED\n");
        printf("--------------------\n");
    }

    /*  Vec4 x1 = {
         floor((r_ij_x) / table_step) + table_size / 2,
         floor((r_ij_y) / table_step) + table_size / 2,
         floor((r_ij_z) / table_step) + table_size / 2};
     Vec4 x2 = {
         floor((r_ij_x) / table_step) + table_size / 2 + 1,
         floor((r_ij_y) / table_step) + table_size / 2,
         floor((r_ij_z) / table_step) + table_size / 2};
     Vec4 x3 = {
         floor((r_ij_x) / table_step) + table_size / 2,
         floor((r_ij_y) / table_step) + table_size / 2 + 1,
         floor((r_ij_z) / table_step) + table_size / 2};
     Vec4 x4 = {
         floor((r_ij_x) / table_step) + table_size / 2 + 1,
         floor((r_ij_y) / table_step) + table_size / 2 + 1,
         floor((r_ij_z) / table_step) + table_size / 2};
     Vec4 x5 = {
         floor((r_ij_x) / table_step) + table_size / 2,
         floor((r_ij_y) / table_step) + table_size / 2,
         floor((r_ij_z) / table_step) + table_size / 2 + 1};
     Vec4 x6 = {
         floor((r_ij_x) / table_step) + table_size / 2 + 1,
         floor((r_ij_y) / table_step) + table_size / 2,
         floor((r_ij_z) / table_step) + table_size / 2 + 1};
     Vec4 x7 = {
         floor((r_ij_x) / table_step) + table_size / 2,
         floor((r_ij_y) / table_step) + table_size / 2 + 1,
         floor((r_ij_z) / table_step) + table_size / 2 + 1};
     Vec4 x8 = {
         floor((r_ij_x) / table_step) + table_size / 2 + 1,
         floor((r_ij_y) / table_step) + table_size / 2 + 1,
         floor((r_ij_z) / table_step) + table_size / 2 + 1};

     Vec3 pt = {
         (r_ij_x) / table_step + table_size / 2,
         (r_ij_y) / table_step + table_size / 2,
         (r_ij_z) / table_step + table_size / 2};

     Vec3 v1 = table[(int)x1.x][(int)x1.y][(int)x1.z];
     Vec3 v2 = table[(int)x2.x][(int)x2.y][(int)x2.z];
     Vec3 v3 = table[(int)x3.x][(int)x3.y][(int)x3.z];
     Vec3 v4 = table[(int)x4.x][(int)x4.y][(int)x4.z];
     Vec3 v5 = table[(int)x5.x][(int)x5.y][(int)x5.z];
     Vec3 v6 = table[(int)x6.x][(int)x6.y][(int)x6.z];
     Vec3 v7 = table[(int)x7.x][(int)x7.y][(int)x7.z];
     Vec3 v8 = table[(int)x8.x][(int)x8.y][(int)x8.z];

     Vec3 result;

     x1.w = v1.x;
     x2.w = v2.x;
     x3.w = v3.x;
     x4.w = v4.x;
     x5.w = v5.x;
     x6.w = v6.x;
     x7.w = v7.x;
     x8.w = v8.x;

     result.x = lerp3D(pt.x, pt.y, pt.z, x1, x2, x3, x4, x5, x6, x7, x8);

     x1.w = v1.z;
     x2.w = v2.z;
     x3.w = v3.z;
     x4.w = v4.z;
     x5.w = v5.z;
     x6.w = v6.z;
     x7.w = v7.z;
     x8.w = v8.z;

     result.y = lerp3D(pt.x, pt.y, pt.z, x1, x2, x3, x4, x5, x6, x7, x8);

     x1.w = v1.z;
     x2.w = v2.z;
     x3.w = v3.z;
     x4.w = v4.z;
     x5.w = v5.z;
     x6.w = v6.z;
     x7.w = v7.z;
     x8.w = v8.z;

     result.z = lerp3D(pt.x, pt.y, pt.z, x1, x2, x3, x4, x5, x6, x7, x8);

     return result; */

    return table[(int)floor((r_ij_x) / table_step) + table_size / 2]
                [(int)floor((r_ij_y) / table_step) + table_size / 2]
                [(int)floor((r_ij_z) / table_step) + table_size / 2];
}

Vec3 *edwald_summation(System *system, double *args)
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

    // Bring back all the particles in the (0,0,0) cell
    restore_positions_in_lattice_first_cell(particles, n_particles);

    const int SELECTED_PARTICLE = 0;
    static int r_k_print_countdown = 1;
    double real_sum = 0;
    double k_sum = 0;

    for (size_t i = 0; i < n_particles - 1; i++)
    {
        for (size_t j = i + 1; j < n_particles; j++)
        {

            Vec3 real_space_force = compute_real_space_force(
                particles[i].x - particles[j].x,
                particles[i].y - particles[j].y,
                particles[i].z - particles[j].z);

            /* Vec3 reciprocal_space_force = compute_reciprocal_space_force(
                particles[i].x - particles[j].x,
                particles[i].y - particles[j].y,
                particles[i].z - particles[j].z); */

            Vec3 reciprocal_space_force = tabulated_reciprocal_space_term(
                particles[i].x - particles[j].x,
                particles[i].y - particles[j].y,
                particles[i].z - particles[j].z);

            tmp_forces[i].x += FORCE_TYPE_CONSTANT * particles[i].charge * particles[j].charge * (real_space_force.x + reciprocal_space_force.x);
            tmp_forces[i].y += FORCE_TYPE_CONSTANT * particles[i].charge * particles[j].charge * (real_space_force.y + reciprocal_space_force.y);
            tmp_forces[i].z += FORCE_TYPE_CONSTANT * particles[i].charge * particles[j].charge * (real_space_force.z + reciprocal_space_force.z);

            tmp_forces[j].x -= FORCE_TYPE_CONSTANT * particles[i].charge * particles[j].charge * (real_space_force.x + reciprocal_space_force.x);
            tmp_forces[j].y -= FORCE_TYPE_CONSTANT * particles[i].charge * particles[j].charge * (real_space_force.y + reciprocal_space_force.y);
            tmp_forces[j].z -= FORCE_TYPE_CONSTANT * particles[i].charge * particles[j].charge * (real_space_force.z + reciprocal_space_force.z);

            if (i == SELECTED_PARTICLE)
            {
                k_sum += reciprocal_space_force.x;
                k_sum += reciprocal_space_force.y;
                k_sum += reciprocal_space_force.z;

                real_sum += real_space_force.x;
                real_sum += real_space_force.y;
                real_sum += real_space_force.z;
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