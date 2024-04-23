#include <math.h>
#include <stdlib.h>

#include "constants.h"

int restore_lattice_first_cell(double* pos, int n_particles) {
    for (size_t i = 0; i < n_particles * 3; i += 3) {
        pos[i + 0] = pos[i + 0] - rint(pos[i + 0] / L) * L;
        pos[i + 1] = pos[i + 1] - rint(pos[i + 1] / L) * L;
        pos[i + 2] = pos[i + 2] - rint(pos[i + 2] / L) * L;
    }
    return 0;
}

double* edwald_summation(double t, double* pos, double* vel, int n_particles, double* args, int n_args) {
    double* forces = (double*)calloc(n_particles * 3, sizeof(double));
    if (!forces) {
        return NULL;
    }
    double* f_vec = (double*)calloc(3, sizeof(double));
    if (!f_vec) {
        free(forces);
        return NULL;
    }

    double q1 = 0.1;
    double q2 = 0.1;
    double V = pow(L, 3);
    double sigma = 0.0001;
    int k_range = 1;
    int n = 1;

    // Riporta tutte le particelle nella cella (0,0,0)
    restore_lattice_first_cell(pos, n_particles);

    for (size_t i = 0; i < n_particles * 3; i += 3) {
        // TERMINE SPAZIO REALE

        for (size_t j = 0; j < n_particles * 3; j += 3) {
            if (i == j) continue;  // Evita il conto di una particella con se stessa

            f_vec[0] = pos[i + 0] - pos[j + 0] + 0 * L;  // Componente x vettore direzione forza
            f_vec[1] = pos[i + 1] - pos[j + 1] + 0 * L;  // Componente y vettore direzione forza
            f_vec[2] = pos[i + 2] - pos[j + 2] + 0 * L;  // Componente z vettore direzione forza

            double f_vec_mag = sqrt(pow(f_vec[0], 2) + pow(f_vec[1], 2) + pow(f_vec[2], 2));
            double classical_force = q1 * q2 / (4 * PI * E0 * pow(f_vec_mag, 3));
            double edwald_correction = SQ2 / sqrt(sigma * PI) * exp(-pow(f_vec_mag / sigma, 2) / 2) * f_vec_mag + erfc(f_vec_mag / (sigma * SQ2));
            forces[j + 0] -= f_vec[0] * classical_force * edwald_correction;
            forces[j + 1] -= f_vec[1] * classical_force * edwald_correction;
            forces[j + 2] -= f_vec[2] * classical_force * edwald_correction;
        }

        // TERMINE SPAZIO RECIPROCO
        for (int k_x = -k_range; k_x < k_range; k_x++) {
            for (int k_y = -k_range; k_y < k_range; k_y++) {
                for (int k_z = -k_range; k_z < k_range; k_z++) {
                    if (k_x != 0 && k_y != 0 && k_z != 0) {  // Rimuovo la componente k = 0  che va analizzata separatamente
                        for (size_t j = 0; j < n_particles * 3; j += 3) {
                            if (i == j) continue;
                            double k_vec_mag_square = (k_x * k_x + k_y * k_y + k_z * k_z) * 4 * PI * PI / L / L;
                            double k_vec_pos_dot_prod = (k_x * (pos[i + 0] - pos[j + 0]) + k_y * (pos[i + 1] - pos[j + 1]) + k_z * (pos[i + 2] - pos[j + 2])) * 2 * PI / L;
                            double f_vec_mag = q1 * q2 / (V * E0) * exp(-sigma * sigma * k_vec_mag_square / 2) / k_vec_mag_square * sin(k_vec_pos_dot_prod);
                            forces[j + 0] -= k_x * 2 * PI / L * f_vec_mag;
                            forces[j + 1] -= k_y * 2 * PI / L * f_vec_mag;
                            forces[j + 2] -= k_z * 2 * PI / L * f_vec_mag;
                        }
                    }
                }
            }
        }
    }

    free(f_vec);
    return forces;
}
