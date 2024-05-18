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

double ALPHA = (5.8 / CELL_LENGHT);

// #define ALPHA 1e-8
#define SIGMA (1 / (SQR_2 * ALPHA))
#define N_PARTICLES 300
#define USE_TABULATION_EDWALD_RECIPROCAL_SPACE 1
#define FORCE_TYPE_CONSTANT 1 // F = FORCE_TYPE_CONSTANT * (Q1*Q2)/r^2
#define RECIPROCAL_SPACE_TABLE_SIZE 32
#define LAMBDA (CELL_LENGHT / 2) // Yukawa potential factor (0 = coulomb)
#define USE_COULOMB 1
#define USE_PARALLELIZATION 1
#define N_THREADS 5

int N_K_RANGE = 5;

#endif