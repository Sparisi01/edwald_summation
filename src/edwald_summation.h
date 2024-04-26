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

double* edwald_summation(double t, double* pos, double* vel, int n_particles, double* args, int n_args) {
  double* forces = (double*)calloc(n_particles * 3, sizeof(double));
  if (!forces) {
    return NULL;
  }

  double real_mag;
  double reciprocal_mag;

  // Riporta tutte le particelle nella cella (0,0,0)
  restore_lattice_positions_in_first_cell(pos, n_particles);

  for (size_t i = 0; i < (n_particles - 1) * 3; i += 3) {
    // TERMINE SPAZIO REALE
    ///////////////////////////////////////////////////////////////////////////////////////////
    for (size_t j = i + 3; j < n_particles * 3; j += 3) {
      double r_ij_x = pos[i + 0] - pos[j + 0] - rint(pos[i + 0] - pos[j + 0] / CELL_L) * CELL_L;  // Componente x vettore direzione forza
      double r_ij_y = pos[i + 1] - pos[j + 1] - rint(pos[i + 1] - pos[j + 1] / CELL_L) * CELL_L;  // Componente y vettore direzione forza
      double r_ij_z = pos[i + 2] - pos[j + 2] - rint(pos[i + 2] - pos[j + 2] / CELL_L) * CELL_L;  // Componente z vettore direzione forza
      double r_ij = sqrt(r_ij_x * r_ij_x + r_ij_y * r_ij_y + r_ij_z * r_ij_z);

      if (r_ij >= CELL_L / 2) {
        continue;  // CUTOF potenziale
      }

      double classical_force = Q * -Q / (r_ij * r_ij * r_ij);
      double edwald_correction = 2 * ALPHA / SQR_PI * exp(-(ALPHA * r_ij) * (ALPHA * r_ij)) * r_ij + erfc(ALPHA * r_ij);

      forces[i + 0] += FORCE_TYPE_CONSTANT * r_ij_x * classical_force * edwald_correction;
      forces[i + 1] += FORCE_TYPE_CONSTANT * r_ij_y * classical_force * edwald_correction;
      forces[i + 2] += FORCE_TYPE_CONSTANT * r_ij_z * classical_force * edwald_correction;

      forces[j + 0] -= FORCE_TYPE_CONSTANT * r_ij_x * classical_force * edwald_correction;
      forces[j + 1] -= FORCE_TYPE_CONSTANT * r_ij_y * classical_force * edwald_correction;
      forces[j + 2] -= FORCE_TYPE_CONSTANT * r_ij_z * classical_force * edwald_correction;

      real_mag = r_ij * classical_force * edwald_correction;
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

            reciprocal_mag = sqrt(k_mag_sq) * f_vec_mag;
          }
        }
      }
    }
    ///////////////////////////////////////////////////////////////////////////////////////////
  }

  return forces;
}

double* edwald_summation_table(double t, double* pos, double* vel, int n_particles, double* args, int n_args) {
  double* forces = (double*)calloc(n_particles * 3, sizeof(double));
  if (!forces) {
    return NULL;
  }

  // Riporta tutte le particelle nella cella (0,0,0)
  restore_lattice_positions_in_first_cell(pos, n_particles);

  double* erfcTable = args;

  for (size_t i = 0; i < (n_particles - 1) * 3; i += 3) {
    // TERMINE SPAZIO REALE
    ///////////////////////////////////////////////////////////////////////////////////////////
    for (size_t j = i + 3; j < n_particles * 3; j += 3) {
      double r_ij_x = pos[i + 0] - pos[j + 0] - rint(pos[i + 0] - pos[j + 0] / CELL_L) * CELL_L;  // Componente x vettore direzione forza
      double r_ij_y = pos[i + 1] - pos[j + 1] - rint(pos[i + 1] - pos[j + 1] / CELL_L) * CELL_L;  // Componente y vettore direzione forza
      double r_ij_z = pos[i + 2] - pos[j + 2] - rint(pos[i + 2] - pos[j + 2] / CELL_L) * CELL_L;  // Componente z vettore direzione forza
      double r_ij = sqrt(r_ij_x * r_ij_x + r_ij_y * r_ij_y + r_ij_z * r_ij_z);

      if (r_ij >= CELL_L / 2) {
        continue;  // CUTOF potenziale
      }

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
