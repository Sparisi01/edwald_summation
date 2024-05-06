#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "constants.h"
#include "edwald_summation.h"

int main(int argc, char const *argv[])
{

    FILE *file = fopen("../output/error_table_reciprocal.dat", "w");

    double step = 1e-2;
    int n_step = CELL_LENGHT / step;

    // Generate table
    tabuled_reciprocal_space_term(0, 0, 0);

    /* for (size_t i = 0; i < n_step; i++)
        for (size_t j = 0; j < n_step; j++)
            for (size_t k = 0; k < n_step; k++)
            {
                Vec3 table_value = tabuled_reciprocal_space_term(i * step, j * step, k * step);
                Vec3 exact_value = compute_reciprocal_space_force(i * step, j * step, k * step);

                double diff = sqrt(pow(table_value.x - exact_value.x, 2) +
                                   pow(table_value.y - exact_value.y, 2) +
                                   pow(table_value.z - exact_value.z, 2));

                double table_value_mod = sqrt(table_value.x * table_value.x + table_value.y * table_value.y + table_value.z + table_value.z);
                double exact_value_mod = sqrt(exact_value.x * exact_value.x + exact_value.y * exact_value.y + exact_value.z + exact_value.z);

                fprintf(file, "%.5E %.5E %.5E %.5E\n", i * step, diff, table_value_mod, exact_value_mod);
            } */

    for (int k = -n_step; k < n_step; k++)
    {
        Vec3 table_value = tabulated_reciprocal_space_term(1, 1, k * step);
        Vec3 exact_value = compute_reciprocal_space_force(1, 1, k * step);

        double diff = sqrt(pow(table_value.x - exact_value.x, 2) +
                           pow(table_value.y - exact_value.y, 2) +
                           pow(table_value.z - exact_value.z, 2));

        double table_value_mod = sqrt(table_value.x * table_value.x + table_value.y * table_value.y + table_value.z + table_value.z);
        double exact_value_mod = sqrt(exact_value.x * exact_value.x + exact_value.y * exact_value.y + exact_value.z + exact_value.z);

        fprintf(file, "%.5E %.5E %.5E %.5E\n", k * step, diff, table_value_mod, exact_value_mod);
    }

    return 0;
}
