#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "constants.h"
#include "edwald_summation.h"

int main(int argc, char const *argv[])
{

    FILE *file = fopen("../output/error_table_reciprocal_128_lin.dat", "w");

    double n_step = 100;
    double a = 2 * CELL_LENGHT / (RECIPROCAL_SPACE_TABLE_SIZE - 4);
    double step = a / n_step;

    // Generate table
    tabulated_reciprocal_space_term(0, 0, 0);

    for (int k = -RECIPROCAL_SPACE_TABLE_SIZE / 2 * n_step; k < (RECIPROCAL_SPACE_TABLE_SIZE / 2 - 1) * n_step; k++)
    {
        Vec3 table_value = tabulated_reciprocal_space_term(a, a, k * step);
        Vec3 exact_value = compute_reciprocal_space_force(a, a, k * step);

        fprintf(file, "%.5E %.5E %.5E %.5E\n", k * step, table_value.z / exact_value.z, table_value.z, exact_value.z);
    }

    return 0;
}
