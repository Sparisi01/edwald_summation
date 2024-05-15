#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "constants.h"
#include "edwald_summation.h"
#include "statistic.h"

#define N_SAMPLES 10000

int main(int argc, char const *argv[])
{

    srand(RAND_SEED);
    printf("Starting evaluation for %d samples\n", N_SAMPLES);
    printf("Matrix size: %d\nCell lenght: %f\n", RECIPROCAL_SPACE_TABLE_SIZE, CELL_LENGHT);

    // Generate table
    tabulated_reciprocal_space_term((Vec3){0, 0, 0});

    FILE *table_error_file = fopen("../output/table_error_file.dat", "w");

    double *values_x = (double *)malloc(sizeof(double) * N_SAMPLES);
    double *values_y = (double *)malloc(sizeof(double) * N_SAMPLES);
    double *values_z = (double *)malloc(sizeof(double) * N_SAMPLES);
    for (size_t i = 0; i < N_SAMPLES; i++)
    {
        double rnd1 = (2 * rand() / (RAND_MAX + 1.) - 1) * CELL_LENGHT;
        double rnd2 = (2 * rand() / (RAND_MAX + 1.) - 1) * CELL_LENGHT;
        double rnd3 = (2 * rand() / (RAND_MAX + 1.) - 1) * CELL_LENGHT;

        Vec3 table_value = tabulated_reciprocal_space_term((Vec3){rnd1, rnd2, rnd3});
        Vec3 exact_value = compute_reciprocal_space_force((Vec3){rnd1, rnd2, rnd3});

        values_x[i] = (table_value.x - exact_value.x) / exact_value.x;
        values_y[i] = (table_value.y - exact_value.y) / exact_value.y;
        values_z[i] = (table_value.z - exact_value.z) / exact_value.z;

        fprintf(table_error_file, "%lf %lf %lf %lf %lf %lf %lf\n", sqrt(rnd1 * rnd1 + rnd2 * rnd2 + rnd3 * rnd3), values_x[i], values_y[i], values_z[i], exact_value.x, exact_value.y, exact_value.z);
    }

    printf("Relative error x: %.3E ± %.3E\n", mean(values_x, 0, N_SAMPLES - 1), stddev(values_x, 0, N_SAMPLES - 1));
    printf("Relative error y: %.3E ± %.3E\n", mean(values_y, 0, N_SAMPLES - 1), stddev(values_y, 0, N_SAMPLES - 1));
    printf("Relative error z: %.3E ± %.3E\n", mean(values_z, 0, N_SAMPLES - 1), stddev(values_z, 0, N_SAMPLES - 1));

    fclose(table_error_file);
    return 0;
}

// Matrix size: 64
//  Cell lenght: 2.000000
//  STARTING GENERATION TABLE
//  TABLE COMPLETED
//  --------------------
//  Relative error x: -1.624E-03 ± 3.698E-03
//  Relative error y: -1.415E-03 ± 3.779E-03
//  Relative error z: -1.554E-03 ± 3.758E-03
