#ifndef CONSTANTS_H
#define CONSTANTS_H

#define PI 3.1415926535897932384626433 // An over-accurate PI
#define SQR_2 1.4142135623730950488    // Square root of 2
#define SQR_3 1.73205080757            // Square root of 3
#define SQR_PI 1.77245385091           // Square root of PI
#define CELL_LENGHT 2.
#define CELL_VOLUME (CELL_LENGHT * CELL_LENGHT * CELL_LENGHT)
#define RAND_SEED 3
#define SIGMA_VELOCITIES (sqrt(1))
#define ALPHA (2 / CELL_LENGHT)
// #define ALPHA 1e8
#define SIGMA (1 / (SQR_2 * ALPHA))
#define N_PARTICLES 20
#define USE_TABULATION_EDWALD_RECIPROCAL_SPACE 0
#define N_K_RANGE 5
#define FORCE_TYPE_CONSTANT 1 // F = FORCE_TYPE_CONSTANT * (Q1*Q2)/r^2
#define RECIPROCAL_SPACE_TABLE_SIZE 128
#define LAMBDA 0 // Yukawa potential factor (0 = coulomb)
#define USE_PARALLELIZATION 0
#define N_THREADS 10
#endif