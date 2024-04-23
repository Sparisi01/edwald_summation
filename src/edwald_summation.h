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

    double q1 = 10;
    double q2 = 10;
    double V = pow(L, 3);
    double sigma = 1e-1;
    int n_k_range = 0;
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

            forces[i + 0] += f_vec[0] * classical_force * edwald_correction;
            forces[i + 1] += f_vec[1] * classical_force * edwald_correction;
            forces[i + 2] += f_vec[2] * classical_force * edwald_correction;
        }

        // TERMINE SPAZIO RECIPROCO
        for (int n_k_x = -n_k_range; n_k_x < n_k_range; n_k_x++) {
            for (int n_k_y = -n_k_range; n_k_y < n_k_range; n_k_y++) {
                for (int n_k_z = -n_k_range; n_k_z < n_k_range; n_k_z++) {
                    if (n_k_x != 0 && n_k_y != 0 && n_k_z != 0) {
                        // Rimuovo la componente k = 0  che va analizzata separatamente
                        double k_x = n_k_x * 2 * PI / L;  // Componente x vettore spazione reciproco
                        double k_y = n_k_y * 2 * PI / L;  // Componente y vettore spazione reciproco
                        double k_z = n_k_z * 2 * PI / L;  // Componente z vettore spazione reciproco
                        for (size_t j = 0; j < n_particles * 3; j += 3) {
                            if (i == j) continue;
                            double k_vec_mag_square = (k_x * k_x + k_y * k_y + k_z * k_z);
                            double k_vec_pos_dot_prod = (k_x * (pos[i + 0] - pos[j + 0]) + k_y * (pos[i + 1] - pos[j + 1]) + k_z * (pos[i + 2] - pos[j + 2]));
                            double f_vec_mag = q1 * q2 / (V * E0) * exp(-sigma * sigma * k_vec_mag_square / 2) / k_vec_mag_square * sin(k_vec_pos_dot_prod);
                            forces[i + 0] += k_x * f_vec_mag;
                            forces[i + 1] += k_y * f_vec_mag;
                            forces[i + 2] += k_z * f_vec_mag;
                        }
                    }
                }
            }
        }
    }

    free(f_vec);
    return forces;
}
