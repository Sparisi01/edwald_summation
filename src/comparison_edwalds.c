#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "constants.h"
#include "edwald_summation.h"
#include "statistic.h"

#define N_SAMPLES 100000

int main(int argc, char const *argv[])
{

    FILE *file = fopen("../output/error_table_reciprocal_128_lin.dat", "w");

    double n_step = 100;
    double a = 2 * CELL_LENGHT / (RECIPROCAL_SPACE_TABLE_SIZE - 4);
    double step = a / n_step;

    printf("Starting evaluation for %d samples\n", N_SAMPLES);
    printf("Matrix size: %d\nCell lenght: %f\n", RECIPROCAL_SPACE_TABLE_SIZE, CELL_LENGHT);

    // Generate table
    tabulated_reciprocal_space_term(0, 0, 0);

    /* for (int k = -RECIPROCAL_SPACE_TABLE_SIZE / 2 * n_step; k < (RECIPROCAL_SPACE_TABLE_SIZE / 2 - 1) * n_step; k++)
    {
        Vec3 table_value = tabulated_reciprocal_space_term(a, a, k * step);
        Vec3 exact_value = compute_reciprocal_space_force(a, a, k * step);

        fprintf(file, "%.5E %.5E %.5E %.5E\n", k * step, table_value.z / exact_value.z, table_value.z, exact_value.z);
    } */

    double *values_x = (double *)malloc(sizeof(double) * N_SAMPLES);
    double *values_y = (double *)malloc(sizeof(double) * N_SAMPLES);
    double *values_z = (double *)malloc(sizeof(double) * N_SAMPLES);
    for (size_t i = 0; i < N_SAMPLES; i++)
    {
        double rnd1 = (2 * rand() / (RAND_MAX + 1.) - 1) * CELL_LENGHT;
        double rnd2 = (2 * rand() / (RAND_MAX + 1.) - 1) * CELL_LENGHT;
        double rnd3 = (2 * rand() / (RAND_MAX + 1.) - 1) * CELL_LENGHT;
        Vec3 table_value = tabulated_reciprocal_space_term(rnd1, rnd2, rnd3);
        Vec3 exact_value = compute_reciprocal_space_force(rnd1, rnd2, rnd3);

        values_x[i] = (table_value.x - exact_value.x) / exact_value.x;
        values_y[i] = (table_value.y - exact_value.y) / exact_value.y;
        values_z[i] = (table_value.z - exact_value.z) / exact_value.z;
    }

    printf("Relative error x: %.3E ± %.3E\n", mean(values_x, 0, N_SAMPLES - 1), stddev(values_x, 0, N_SAMPLES - 1));
    printf("Relative error y: %.3E ± %.3E\n", mean(values_y, 0, N_SAMPLES - 1), stddev(values_y, 0, N_SAMPLES - 1));
    printf("Relative error z: %.3E ± %.3E\n", mean(values_z, 0, N_SAMPLES - 1), stddev(values_z, 0, N_SAMPLES - 1));

    fclose(file);
    return 0;
}
