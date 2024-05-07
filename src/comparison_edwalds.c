#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "constants.h"
#include "edwald_summation.h"
#include "statistic.h"

int main(int argc, char const *argv[])
{

    FILE *file = fopen("../output/error_table_reciprocal_128_lin.dat", "w");

    double n_step = 100;
    double a = 2 * CELL_LENGHT / (RECIPROCAL_SPACE_TABLE_SIZE - 4);
    double step = a / n_step;

    // Generate table
    tabulated_reciprocal_space_term(0, 0, 0);

    /* for (int k = -RECIPROCAL_SPACE_TABLE_SIZE / 2 * n_step; k < (RECIPROCAL_SPACE_TABLE_SIZE / 2 - 1) * n_step; k++)
    {
        Vec3 table_value = tabulated_reciprocal_space_term(a, a, k * step);
        Vec3 exact_value = compute_reciprocal_space_force(a, a, k * step);

        fprintf(file, "%.5E %.5E %.5E %.5E\n", k * step, table_value.z / exact_value.z, table_value.z, exact_value.z);
    } */

    double values[10000];
    for (size_t i = 0; i < 10000; i++)
    {
        double rnd1 = (rand() / (RAND_MAX + 1.) - 2) * CELL_LENGHT;
        double rnd2 = (rand() / (RAND_MAX + 1.) - 2) * CELL_LENGHT;
        double rnd3 = (rand() / (RAND_MAX + 1.) - 2) * CELL_LENGHT;
        Vec3 table_value = tabulated_reciprocal_space_term(rnd1, rnd2, rnd3);
        Vec3 exact_value = compute_reciprocal_space_force(rnd1, rnd2, rnd3);

        if (exact_value.x != 0)
            values[i] = table_value.x / exact_value.x;
    }

    double m = mean(values, 0, 9999);
    double d = stddev(values, 0, 9999);

    printf("%lf ± %lf", m, d);

    fclose(file);
    return 0;
}
