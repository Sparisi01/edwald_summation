#define PI 3.1415926535897932384626433
#define SQR_2 1.4142135623730950488
#define SQR_3 1.73205080757
#define SQR_PI 1.77245385091
#define CELL_L 2                                // Cell length
#define CELL_V (CELL_L * CELL_L * CELL_L)       // Cell volume
#define N_PARTICLES 270                         //
#define ALPHA (5.6 / CELL_L)                    //
#define SIGMA (1 / (SQR_2 * ALPHA))             // STD deviation of density correction to point particle
#define N_K_RANGE 2                             // Range of reciprocal lattice cell computed
#define FORCE_TYPE_CONSTANT 1                   // F = FORCE_TYPE_CONSTANT * (Q1*Q2)/r^2
#define ERFC_TABLE_IN 0                         // Start value erf table generator
#define ERFC_TABLE_FIN (ALPHA * CELL_L * SQR_3) // End value erf table generator
#define ERFC_TABLE_PRECISION 1E-5               //
#define SEED 5                                  // Seed rand() function, made for replicate the same system N times
#define Q 0.01