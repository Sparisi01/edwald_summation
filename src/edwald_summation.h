#include <math.h>
#include <stdlib.h>

#include "constants.h"

int restore_lattice_positions_in_first_cell(double* pos, int n_particles) {
    for (size_t i = 0; i < n_particles * 3; i += 3) {
        pos[i + 0] = pos[i + 0] - rint(pos[i + 0] / CELL_L) * CELL_L;
        pos[i + 1] = pos[i + 1] - rint(pos[i + 1] / CELL_L) * CELL_L;
        pos[i + 2] = pos[i + 2] - rint(pos[i + 2] / CELL_L) * CELL_L;
    }
    return 0;
}

double* edwald_summation(double t, double* pos, double* vel, int n_particles, void* args) {
    double* forces = (double*)calloc(n_particles * 3, sizeof(double));
    if (!forces) {
        return NULL;
    }

    // Array containing particles charges
    double* charge_particle_array = (double*)args;
    // Bring back all the particles in the (0,0,0) cell
    restore_lattice_positions_in_first_cell(pos, n_particles);

    for (size_t i = 0; i < (n_particles - 1) * 3; i += 3) {
        // SECTION - Real space force
        for (size_t j = i + 3; j < n_particles * 3; j += 3) {
            /**NOTE - all the following variable are coordinate system independet.
             * NOTE - se applico il cutof del potenziale, fare un ciclo su le celle nello spazio
             * reale oppure applicare direttamente r_int per trovare la coppia più vicina è equivalente.
             **/
            double r_ij_x = pos[i + 0] - pos[j + 0] - rint((pos[i + 0] - pos[j + 0]) / CELL_L) * CELL_L;  //
            double r_ij_y = pos[i + 1] - pos[j + 1] - rint((pos[i + 1] - pos[j + 1]) / CELL_L) * CELL_L;  //
            double r_ij_z = pos[i + 2] - pos[j + 2] - rint((pos[i + 2] - pos[j + 2]) / CELL_L) * CELL_L;  //
            // double r_ij_x = pos[i + 0] - pos[j + 0];  //
            // double r_ij_y = pos[i + 1] - pos[j + 1];  //
            // double r_ij_z = pos[i + 2] - pos[j + 2];  //

            double r_ij = sqrt(r_ij_x * r_ij_x + r_ij_y * r_ij_y + r_ij_z * r_ij_z);  // |r_i - r_j|

            /* if (r_ij >= CELL_L / 2) {
                continue;  // CUTOF potential
            } */

            double classical_force = charge_particle_array[i / 3] * charge_particle_array[j / 3] / (r_ij * r_ij * r_ij);        // Classic Coulomb Like force
            double edwald_correction = 2 * ALPHA / SQR_PI * exp(-(ALPHA * r_ij) * (ALPHA * r_ij)) * r_ij + erfc(ALPHA * r_ij);  // Edwald correction in real space

            // Force on particle i due to j
            forces[i + 0] += FORCE_TYPE_CONSTANT * r_ij_x * classical_force * edwald_correction;
            forces[i + 1] += FORCE_TYPE_CONSTANT * r_ij_y * classical_force * edwald_correction;
            forces[i + 2] += FORCE_TYPE_CONSTANT * r_ij_z * classical_force * edwald_correction;
            // Forces on particle j due to i using third principle
            forces[j + 0] -= FORCE_TYPE_CONSTANT * r_ij_x * classical_force * edwald_correction;
            forces[j + 1] -= FORCE_TYPE_CONSTANT * r_ij_y * classical_force * edwald_correction;
            forces[j + 2] -= FORCE_TYPE_CONSTANT * r_ij_z * classical_force * edwald_correction;
        }

        /*!SECTION */

        // SECTION - Reciprocal space force
        for (int n_k_x = -N_K_RANGE; n_k_x <= N_K_RANGE; n_k_x++) {
            for (int n_k_y = -N_K_RANGE; n_k_y <= N_K_RANGE; n_k_y++) {
                for (int n_k_z = -N_K_RANGE; n_k_z <= N_K_RANGE; n_k_z++) {
                    if (n_k_x == 0 && n_k_y == 0 && n_k_z == 0) {
                        // TODO:
                        continue;  // Rimuovo la componente (0,0,0) che va analizzata separatamente
                    }

                    double k_x = n_k_x * 2 * PI / CELL_L;  // Componente x vettore spazione reciproco
                    double k_y = n_k_y * 2 * PI / CELL_L;  // Componente y vettore spazione reciproco
                    double k_z = n_k_z * 2 * PI / CELL_L;  // Componente z vettore spazione reciproco

                    /**NOTE - La formula non prevede l'esclusione dell'interazione tra particella i e j
                     * se esse appartengono a celle diverse. Noto però che in questo caso il sin da valore
                     * nullo. Posso escludere l'interazione.
                     **/
                    for (size_t j = i + 3; j < n_particles * 3; j += 3) {
                        double k_mag_sq = (k_x * k_x + k_y * k_y + k_z * k_z);  // Magnitude square of reciprocal-lattice vector
                        double dot_prod_pos_k = (k_x * (pos[i + 0] - pos[j + 0]) + k_y * (pos[i + 1] - pos[j + 1]) + k_z * (pos[i + 2] - pos[j + 2]));
                        double f_vec_mag = charge_particle_array[i / 3] * charge_particle_array[j / 3] * 4 * PI / CELL_V * exp(-k_mag_sq / (2 * ALPHA * ALPHA)) / k_mag_sq * sin(dot_prod_pos_k);
                        // Forces on particle i due to j copies
                        forces[i + 0] += FORCE_TYPE_CONSTANT * k_x * f_vec_mag;
                        forces[i + 1] += FORCE_TYPE_CONSTANT * k_y * f_vec_mag;
                        forces[i + 2] += FORCE_TYPE_CONSTANT * k_z * f_vec_mag;
                        // Forces on particle j due to i copies using third principle
                        forces[j + 0] -= FORCE_TYPE_CONSTANT * k_x * f_vec_mag;
                        forces[j + 1] -= FORCE_TYPE_CONSTANT * k_y * f_vec_mag;
                        forces[j + 2] -= FORCE_TYPE_CONSTANT * k_z * f_vec_mag;
                    }
                }
            }
        }
        /*!SECTION */
    }

    return forces;
}

double* edwald_summation_table(double t, double* pos, double* vel, int n_particles, void* args) {
    double* forces = (double*)calloc(n_particles * 3, sizeof(double));
    if (!forces) {
        return NULL;
    }

    // Riporta tutte le particelle nella cella (0,0,0)
    restore_lattice_positions_in_first_cell(pos, n_particles);

    double* erfcTable = (double*)args;

    for (size_t i = 0; i < (n_particles - 1) * 3; i += 3) {
        // TERMINE SPAZIO REALE
        ///////////////////////////////////////////////////////////////////////////////////////////
        for (size_t j = i + 3; j < n_particles * 3; j += 3) {
            double r_ij_x = pos[i + 0] - pos[j + 0] - rint(pos[i + 0] - pos[j + 0] / CELL_L) * CELL_L;  // Componente x vettore direzione forza
            double r_ij_y = pos[i + 1] - pos[j + 1] - rint(pos[i + 1] - pos[j + 1] / CELL_L) * CELL_L;  // Componente y vettore direzione forza
            double r_ij_z = pos[i + 2] - pos[j + 2] - rint(pos[i + 2] - pos[j + 2] / CELL_L) * CELL_L;  // Componente z vettore direzione forza
            double r_ij = sqrt(r_ij_x * r_ij_x + r_ij_y * r_ij_y + r_ij_z * r_ij_z);

            /* if (r_ij >= CELL_L / 2) {
                continue;  // CUTOF potenziale
            } */

            double curr_erfc = (erfcTable[(int)floor(ALPHA * r_ij / ERFC_TABLE_PRECISION) + 1] + erfcTable[(int)floor(ALPHA * r_ij / ERFC_TABLE_PRECISION)]) / 2;

            double classical_force = Q * -Q / (r_ij * r_ij * r_ij);
            double edwald_correction = 2 * ALPHA / SQR_PI * exp(-(ALPHA * r_ij) * (ALPHA * r_ij)) * r_ij + curr_erfc;

            forces[i + 0] += FORCE_TYPE_CONSTANT * r_ij_x * classical_force * edwald_correction;
            forces[i + 1] += FORCE_TYPE_CONSTANT * r_ij_y * classical_force * edwald_correction;
            forces[i + 2] += FORCE_TYPE_CONSTANT * r_ij_z * classical_force * edwald_correction;

            forces[j + 0] -= FORCE_TYPE_CONSTANT * r_ij_x * classical_force * edwald_correction;
            forces[j + 1] -= FORCE_TYPE_CONSTANT * r_ij_y * classical_force * edwald_correction;
            forces[j + 2] -= FORCE_TYPE_CONSTANT * r_ij_z * classical_force * edwald_correction;
        }

        // TERMINE SPAZIO RECIPROCO
        ///////////////////////////////////////////////////////////////////////////////////////////
        for (int n_k_x = -N_K_RANGE; n_k_x <= N_K_RANGE; n_k_x++) {
            for (int n_k_y = -N_K_RANGE; n_k_y <= N_K_RANGE; n_k_y++) {
                for (int n_k_z = -N_K_RANGE; n_k_z <= N_K_RANGE; n_k_z++) {
                    if (n_k_x == 0 && n_k_y == 0 && n_k_z == 0) {
                        // TODO: implementare spazio reciproco cella 0
                        continue;  // Rimuovo la componente (0,0,0) che va analizzata separatamente
                    }
                    double k_x = n_k_x * 2 * PI / CELL_L;  // Componente x vettore spazione reciproco
                    double k_y = n_k_y * 2 * PI / CELL_L;  // Componente y vettore spazione reciproco
                    double k_z = n_k_z * 2 * PI / CELL_L;  // Componente z vettore spazione reciproco
                    for (size_t j = i + 3; j < n_particles * 3; j += 3) {
                        double k_mag_sq = (k_x * k_x + k_y * k_y + k_z * k_z);
                        double dot_prod_pos_k = (k_x * (pos[i + 0] - pos[j + 0]) + k_y * (pos[i + 1] - pos[j + 1]) + k_z * (pos[i + 2] - pos[j + 2]));
                        double f_vec_mag = Q * -Q * 4 * PI / CELL_V * exp(-k_mag_sq / (ALPHA * ALPHA)) / k_mag_sq * sin(dot_prod_pos_k);
                        forces[i + 0] += FORCE_TYPE_CONSTANT * k_y * f_vec_mag;
                        forces[i + 1] += FORCE_TYPE_CONSTANT * k_y * f_vec_mag;
                        forces[i + 2] += FORCE_TYPE_CONSTANT * k_z * f_vec_mag;

                        forces[j + 0] -= FORCE_TYPE_CONSTANT * k_x * f_vec_mag;
                        forces[j + 1] -= FORCE_TYPE_CONSTANT * k_y * f_vec_mag;
                        forces[j + 2] -= FORCE_TYPE_CONSTANT * k_z * f_vec_mag;
                    }
                }
            }
        }
        ///////////////////////////////////////////////////////////////////////////////////////////
    }

    return forces;
}