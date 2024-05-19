#ifndef EDWALD_RECIPROCAL_SPACE_H
#define EDWALD_RECIPROCAL_SPACE_H

#include "../complex_utils.h"
#include "../coulomb/coulomb_potential.h"
#include "../structures.h"
#include <complex.h>
#include <math.h>
#include <stdlib.h>

int _K_RANGE = 5;

double reciprocal_space_potential_coulomb(System *s, double ALPHA)
{

    double base_frequency = (2 * PI / s->cell_lenght);

    //____________________________________________
    // COMPUTE STRUCTURAL FACTOR FOR ALL K-VECTORS

    double complex structure_factor[2 * _K_RANGE + 1][2 * _K_RANGE + 1][2 * _K_RANGE + 1];

    for (int k_x = -_K_RANGE; k_x <= _K_RANGE; k_x++)
    {
        for (int k_y = -_K_RANGE; k_y <= _K_RANGE; k_y++)
        {
            for (int k_z = -_K_RANGE; k_z <= _K_RANGE; k_z++)
            {
                structure_factor[k_x + _K_RANGE][k_y + _K_RANGE][k_z + _K_RANGE] = 0 + 0 * I;
            }
        }
    }

    for (size_t i = 0; i < s->n_particles; i++)
    {
        // Compute exp(i*2π/L*rᵢ) for x,y and z
        double complex e_x = cexp(I * base_frequency * s->particles[i].x);
        double complex e_y = cexp(I * base_frequency * s->particles[i].y);
        double complex e_z = cexp(I * base_frequency * s->particles[i].z);

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
                    structure_factor[k_x + _K_RANGE][k_y + _K_RANGE][k_z + _K_RANGE] +=
                        precompure_pow_x[k_x + _K_RANGE] *
                        precompure_pow_y[k_y + _K_RANGE] *
                        precompure_pow_z[k_z + _K_RANGE] * s->particles[i].charge;
                }
            }
        }
    }

    //________________________________
    // COMPUTE RECIPROCAL SPACE ENERGY

    double sum = 0;
    double volume = pow(s->cell_lenght, 3);

    for (int k_x = -_K_RANGE; k_x <= _K_RANGE; k_x++)
    {
        for (int k_y = -_K_RANGE; k_y <= _K_RANGE; k_y++)
        {
            for (int k_z = -_K_RANGE; k_z <= _K_RANGE; k_z++)
            {
                // If the system is charge neutral ignore the first cell
                if (k_x == 0 && k_y == 0 && k_z == 0) continue;

                double k_mod2 = (k_x * k_x + k_y * k_y + k_z * k_z) * base_frequency * base_frequency;
                double str_quad = (structure_factor[k_x + _K_RANGE][k_y + _K_RANGE][k_z + _K_RANGE] * conj(structure_factor[k_x + _K_RANGE][k_y + _K_RANGE][k_z + _K_RANGE]));
                sum += exp(-(k_mod2 / 4 * ALPHA * ALPHA)) / k_mod2 * str_quad;
            }
        }
    }

    return (4 * PI / volume) * 0.5 * sum;
}

Vec3 reciprocal_space_force()
{
}

#endif