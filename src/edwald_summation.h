#include <math.h>
#include <stdlib.h>

#define PI 3.14159265359
#define SQ2 1.41421356237

double* edwald_summation(double t, double* pos, double* vel, int n_particles, double* args, int n_args) {
    double* forces = (double*)calloc(n_particles * 3, sizeof(double));
    if (!forces) {
        return NULL;
    }

    double* f_vec = (double*)calloc(3, sizeof(double));
    double q1 = 1;
    double q2 = 2;
    double L = 1;
    double V = pow(L, 3);
    double sigma = 1;
    double E0 = 1;
    int i = 0;

    // TERMINE SPAZIO REALE
    int n = 1;
    for (size_t j = 0; j < n_particles * 3; j += 3) {
        f_vec[0] = pos[i + 0] - pos[j + 0] + 0 * L;  // Componente x vettore direzione forza
        f_vec[1] = pos[i + 1] - pos[j + 1] + 0 * L;  // Componente y vettore direzione forza
        f_vec[2] = pos[i + 2] - pos[j + 2] + 0 * L;  // Componente z vettore direzione forza

        double f_vec_mag = sqrt(pow(f_vec[0], 2) + pow(f_vec[1], 2) + pow(f_vec[2], 2));
        double classical_force = q1 * q2 / (4 * PI * E0 * pow(f_vec_mag, 3));
        double edwald_correction = SQ2 / sqrt(sigma * PI) * exp(-pow(f_vec_mag / sigma, 2) / 2) * f_vec_mag + erfc(f_vec_mag / (sigma * SQ2));
        forces[j + 0] += f_vec[0] * classical_force * edwald_correction;
        forces[j + 1] += f_vec[1] * classical_force * edwald_correction;
        forces[j + 2] += f_vec[2] * classical_force * edwald_correction;
    }

    // TERMINE SPAZIO RECIPROCO
    int k_range = 2;

    for (size_t k_x = 0; k_x < k_range; k_x++) {
        for (size_t k_y = 0; k_y < k_range; k_y++) {
            for (size_t k_z = 0; k_z < k_range; k_z++) {
                
            }
        }
    }

    free(f_vec);
    return forces;
}
