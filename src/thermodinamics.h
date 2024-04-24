#include <math.h>

#include "constants.h"

double temperature(double* vel) {
  return 1;
}

double kinetic_energy(double* vel, double* masses, int n_particles) {
  double kinetic_energy = 0;

  for (size_t i = 0; i < n_particles * 3; i += 3) {
    kinetic_energy += 0.5 * masses[i / 3] * (vel[i + 0] * vel[i + 0] + vel[i + 1] * vel[i + 1] + vel[i + 2] * vel[i + 2]);
  }

  return kinetic_energy;
}

double coulomb_potential_energy(double* pos, int n_particles) {
  double q1 = 0.5;
  double q2 = 0.5;
  double pot_energy = 0;

  for (size_t i = 0; i < n_particles * 3; i += 3) {
    for (size_t j = i + 3; j < n_particles * 3; j += 3) {
      double vec_dist_mag = sqrt(pow(pos[i + 0] - pos[j + 0], 2) + pow(pos[i + 1] - pos[j + 1], 2) + pow(pos[i + 2] - pos[j + 2], 2));
      pot_energy += q1 * q2 / (4 * PI * E0 * vec_dist_mag);
    }
  }
  return pot_energy;
}
