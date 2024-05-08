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

/// @brief Function that use a 3D pre-computed vectorial field to perform a
/// trilinear interpolation and estimate the edwald force in r_ij.
/// The first time the function is called the table is filled/loaded from file.
///
/// UPPERBOUND ERROR FOR DIFFERENT MATRIX SIZES
/// 32x32x32: 2.5%
/// 64x64x64: 0.6%
/// 128x128x128: 0.15%
///
/// @param r_ij_x
/// @param r_ij_y
/// @param r_ij_z
/// @return Vec3 containing the estimated force-vector in r_ij
Vec3 tabulated_reciprocal_space_term(double r_ij_x, double r_ij_y, double r_ij_z)
{
    const int table_size = RECIPROCAL_SPACE_TABLE_SIZE;
    // -4 done to have a layer of 4 cells unreachable by the system
    // but used in case of a tricubic interpolation
    // table range is [-(2*CELL_LENGHT + 2*table_step), (2*CELL_LENGHT + 2*table_step)]
    const double table_step = 2 * CELL_LENGHT / (table_size - 4.);
    static int has_been_tabled = 0;
    /* static Vec3 table
        [RECIPROCAL_SPACE_TABLE_SIZE]
        [RECIPROCAL_SPACE_TABLE_SIZE]
        [RECIPROCAL_SPACE_TABLE_SIZE]; */

    static Vec3 ***table;

    if (!has_been_tabled)
    {
        if (table_size % 2 == 1)
        {
            perror("RECIPROCAL_SPACE_TABLE_SIZE must be even");
            exit(EXIT_FAILURE);
        }

        // Alloc 3D vector field
        table = (Vec3 ***)malloc(sizeof(Vec3 **) * table_size);
        if (!table)
        {
            perror("Error allocating table");
            exit(EXIT_FAILURE);
        }
        for (size_t i = 0; i < table_size; i++)
        {
            table[i] = (Vec3 **)malloc(sizeof(Vec3 *) * table_size);
            if (!table[i])
            {
                perror("Error allocating table");
                exit(EXIT_FAILURE);
            }
            for (size_t j = 0; j < table_size; j++)
            {
                table[i][j] = (Vec3 *)malloc(sizeof(Vec3) * table_size);
                if (!table[i][j])
                {
                    perror("Error allocating table");
                    exit(EXIT_FAILURE);
                }
            }
        }

        printf("STARTING GENERATION TABLE\n");
        for (int i = 0; i < table_size; i++)
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

        if (0)
        {
            printf("SAVING TABLE\n");
            FILE *table_file = fopen("./tables/table_L2_N64.dat", "w");
            for (int i = 0; i < table_size; i++)
                for (int j = 0; j < table_size; j++)
                    for (int k = 0; k < table_size; k++)
                        fprintf(table_file,
                                "%d %d %d %.6E %.6E %.6E\n", i, j, k,
                                table[i][j][k].x,
                                table[i][j][k].y,
                                table[i][j][k].z);
            printf("--------------------\n");
        }
    }

    // The particle is contained in a table cell, the index of the box bottom vertex
    // is x0,y0,z0. The value of edwald can be obtained for that vertex by retrieving
    // table[x0][y0][z0]
    int x0 = floor((r_ij_x) / table_step) + table_size / 2;
    int y0 = floor((r_ij_y) / table_step) + table_size / 2;
    int z0 = floor((r_ij_z) / table_step) + table_size / 2;

    // Al coordinate are translated and scaled in order to have a cube
    // of size 1x1x1, the vertex (x0,y0,z0) goes to (0,0,0).
    Vec4 x1 = {0, 0, 0, 0};
    Vec4 x2 = {1, 0, 0, 0};
    Vec4 x3 = {0, 1, 0, 0};
    Vec4 x4 = {1, 1, 0, 0};
    Vec4 x5 = {0, 0, 1, 0};
    Vec4 x6 = {1, 0, 1, 0};
    Vec4 x7 = {0, 1, 1, 0};
    Vec4 x8 = {1, 1, 1, 0};

    // Coordinate of the point in the new coordinate.
    // pt.x ∈ [0,1], pt.y ∈ [0,1], pt.z ∈ [0,1].
    Vec3 pt = {
        (r_ij_x / table_step) - x0 + table_size / 2,
        (r_ij_y / table_step) - y0 + table_size / 2,
        (r_ij_z / table_step) - z0 + table_size / 2,
    };

    // Value of edwald in the vertecies
    Vec3 v1 = table[x0 + 0][y0 + 0][z0 + 0];
    Vec3 v2 = table[x0 + 1][y0 + 0][z0 + 0];
    Vec3 v3 = table[x0 + 0][y0 + 1][z0 + 0];
    Vec3 v4 = table[x0 + 1][y0 + 1][z0 + 0];
    Vec3 v5 = table[x0 + 0][y0 + 0][z0 + 1];
    Vec3 v6 = table[x0 + 1][y0 + 0][z0 + 1];
    Vec3 v7 = table[x0 + 0][y0 + 1][z0 + 1];
    Vec3 v8 = table[x0 + 1][y0 + 1][z0 + 1];

    // Edwald therm is a 3D vectorial field, so we need three interpolations,
    // one for each coordinate.
    Vec3 estimated_force;

    // X FORCE COMPONENT
    x1.w = v1.x;
    x2.w = v2.x;
    x3.w = v3.x;
    x4.w = v4.x;
    x5.w = v5.x;
    x6.w = v6.x;
    x7.w = v7.x;
    x8.w = v8.x;

    estimated_force.x = lerp3D(pt.x, pt.y, pt.z, x1, x2, x3, x4, x5, x6, x7, x8);

    // Y FORCE COMPONENT
    x1.w = v1.y;
    x2.w = v2.y;
    x3.w = v3.y;
    x4.w = v4.y;
    x5.w = v5.y;
    x6.w = v6.y;
    x7.w = v7.y;
    x8.w = v8.y;

    estimated_force.y = lerp3D(pt.x, pt.y, pt.z, x1, x2, x3, x4, x5, x6, x7, x8);

    // Z FORCE COMPONENT
    x1.w = v1.z;
    x2.w = v2.z;
    x3.w = v3.z;
    x4.w = v4.z;
    x5.w = v5.z;
    x6.w = v6.z;
    x7.w = v7.z;
    x8.w = v8.z;

    estimated_force.z = lerp3D(pt.x, pt.y, pt.z, x1, x2, x3, x4, x5, x6, x7, x8);

    return estimated_force;

    /* return table[(int)floor((r_ij_x) / table_step) + table_size / 2]
                [(int)floor((r_ij_y) / table_step) + table_size / 2]
                [(int)floor((r_ij_z) / table_step) + table_size / 2]; */
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

            // Coorrection obtained by the study of mean and stdev
            // on a dataset of trilinear approximation (L=2, table_size = 64).
            // Relative error is:
            // (reciprocal_space_force.a - real_space_force.a) / real_space_force.a
            // with a ∈ {x,y,z}

            // NO CORRECTION:
            // Relative error x: -2.670E-03 ± 3.024E-03
            // Relative error y: -2.643E-03 ± 3.001E-03
            // Relative error z: -2.619E-03 ± 3.079E-03

            // WITH CORRECTION:
            // Relative error x: -4.148E-06 ± 3.032E-03
            // Relative error y : -6.957E-06 ± 3.008E-03
            // Relative error z : -1.372E-05 ± 3.087E-03

            // Code in comparison_edwalds.c

            reciprocal_space_force.x *= 1.00267;
            reciprocal_space_force.y *= 1.00264;
            reciprocal_space_force.z *= 1.00261;

            tmp_forces[i].x += FORCE_TYPE_CONSTANT * particles[i].charge * particles[j].charge * (real_space_force.x + reciprocal_space_force.x);
            tmp_forces[i].y += FORCE_TYPE_CONSTANT * particles[i].charge * particles[j].charge * (real_space_force.y + reciprocal_space_force.y);
            tmp_forces[i].z += FORCE_TYPE_CONSTANT * particles[i].charge * particles[j].charge * (real_space_force.z + reciprocal_space_force.z);

            // Use Newton third principle
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
    }
    if (r_k_print_countdown > 0)
    {
        printf("REAL SUM: %.10E\nRECIPROCAL SUM: %.10E\n", real_sum, k_sum);
        r_k_print_countdown--;
    }

    return tmp_forces;
}
#endif