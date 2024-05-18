#ifndef EDWALD_RECIPROCAL_SPACE_H
#define EDWALD_RECIPROCAL_SPACE_H

#include "../complex_utils.h"
#include "../structures.h"
#include <complex.h>
#include <math.h>
#include <stdlib.h>

int _K_RANGE = 5;

/* double complex ***compute_structure_factor_matrix(System *s)
{
    double complex structure_factor[2 * _K_RANGE + 1][2 * _K_RANGE + 1][2 * _K_RANGE + 1];

    for (size_t i = 0; i < s->n_particles; i++)
    {
        double complex e_x = cos(s->particles[i].x) + I * sin(s->particles[i].x);
        double complex e_y = cos(s->particles[i].y) + I * sin(s->particles[i].y);
        double complex e_z = cos(s->particles[i].z) + I * sin(s->particles[i].z);

        double precompure_pow_x[2 * _K_RANGE + 1];
        double precompure_pow_y[2 * _K_RANGE + 1];
        double precompure_pow_z[2 * _K_RANGE + 1];

        for (int j = -_K_RANGE; j <= _K_RANGE; j++)
        {
            precompure_pow_x[j + _K_RANGE] = cpow(e_x, j);
            precompure_pow_y[j + _K_RANGE] = cpow(e_y, j);
            precompure_pow_z[j + _K_RANGE] = cpow(e_z, j);
        }

        for (int k_x = -_K_RANGE; k_x <= _K_RANGE; k_x++)
        {
            for (int k_y = -_K_RANGE; k_y <= _K_RANGE; k_y++)
            {
                for (int k_z = -_K_RANGE; k_z <= _K_RANGE; k_z++)
                {
                    // Must be removed only if the charge sum is 0
                    // if (k_x == 0 && k_y == 0 && k_z == 0) continue;

                    structure_factor[k_x + _K_RANGE][k_y + _K_RANGE][k_z + _K_RANGE] +=
                        precompure_pow_x[k_x + _K_RANGE] *
                        precompure_pow_y[k_y + _K_RANGE] *
                        precompure_pow_z[k_z + _K_RANGE] * s->particles[i].charge;
                }
            }
        }
    }

    return structure_factor;
} */

double reciprocal_space_potential(System *s, double ALPHA)
{
    double sum = 0;
    double C = 0;

    double complex structure_factor[2 * _K_RANGE + 1][2 * _K_RANGE + 1][2 * _K_RANGE + 1];

    for (size_t i = 0; i < s->n_particles; i++)
    {
        double complex e_x = cos(s->particles[i].x) + I * sin(s->particles[i].x);
        double complex e_y = cos(s->particles[i].y) + I * sin(s->particles[i].y);
        double complex e_z = cos(s->particles[i].z) + I * sin(s->particles[i].z);

        double precompure_pow_x[2 * _K_RANGE + 1];
        double precompure_pow_y[2 * _K_RANGE + 1];
        double precompure_pow_z[2 * _K_RANGE + 1];

        for (int j = -_K_RANGE; j <= _K_RANGE; j++)
        {
            precompure_pow_x[j + _K_RANGE] = cpow(e_x, j);
            precompure_pow_y[j + _K_RANGE] = cpow(e_y, j);
            precompure_pow_z[j + _K_RANGE] = cpow(e_z, j);
        }

        for (int k_x = -_K_RANGE; k_x <= _K_RANGE; k_x++)
        {
            for (int k_y = -_K_RANGE; k_y <= _K_RANGE; k_y++)
            {
                for (int k_z = -_K_RANGE; k_z <= _K_RANGE; k_z++)
                {
                    // Must be removed only if the charge sum is 0
                    // if (k_x == 0 && k_y == 0 && k_z == 0) continue;

                    structure_factor[k_x + _K_RANGE][k_y + _K_RANGE][k_z + _K_RANGE] +=
                        precompure_pow_x[k_x + _K_RANGE] *
                        precompure_pow_y[k_y + _K_RANGE] *
                        precompure_pow_z[k_z + _K_RANGE] * s->particles[i].charge;
                }
            }
        }
    }

    // double complex ***structure_factor = compute_structure_factor_matrix(s);

    for (int k_x = -_K_RANGE; k_x <= _K_RANGE; k_x++)
    {
        for (int k_y = -_K_RANGE; k_y <= _K_RANGE; k_y++)
        {
            for (int k_z = -_K_RANGE; k_z <= _K_RANGE; k_z++)
            {
                if (k_x == 0 && k_y == 0 && k_z == 0) continue;

                double k_mod2 = k_x * k_x + k_y * k_y + k_z * k_z;
                double tmp_exp = exp(-k_mod2 / (4 * ALPHA)) / k_mod2;
                sum += tmp_exp * pow(cabs(structure_factor[k_x + _K_RANGE][k_y + _K_RANGE][k_z + _K_RANGE]), 2);
                C += tmp_exp;
            }
        }
    }

    return (4 * PI / pow(s->cell_lenght, 3)) * 0.5 * (sum - C);
}

Vec3 reciprocal_space_force()
{
}

#endif