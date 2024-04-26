#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "constants.h"
#include "edwald_summation.h"
#include "thermodinamics.h"
#include "verlet_propagation.h"

// using namespace std;

/**
 * TODO: main todo
 * [] RICERCA SIGMA SENSATO
 * [] VERIFICA SISTEMI NOTI
 * [] CUTOF POTENZIALE
 */

int generateErfcTable(double in, double fin, double precision, FILE *tableFile) {
  if (!tableFile) {
    return 1;
  }

  double N = (fin - in) / precision;
  for (size_t i = 0; i < N; i++) {
    fprintf(tableFile, "%lf", erfc(in + i * precision));
  }

  return 0;
}

int loadErfcTable(double in, double fin, double interval, FILE *tableFile) {
  return 1;
}

int main(int argc, char const *argv[]) {
  srand(1);

  double *pos_array = (double *)calloc(N_PARTICLES * 3, sizeof(double));
  double *vel_array = (double *)calloc(N_PARTICLES * 3, sizeof(double));
  double **forces_array_ptr = (double **)malloc(sizeof(double *));
  *forces_array_ptr = (double *)calloc(N_PARTICLES * 3, sizeof(double));
  double *masses_array = (double *)calloc(N_PARTICLES, sizeof(double));

  if (!pos_array || !vel_array || !forces_array_ptr || !*forces_array_ptr || !masses_array) {
    printf("Errore inizializzazione array\n");
    return 1;
  }

  for (size_t i = 0; i < N_PARTICLES; i++) {
    masses_array[i] = 1;
  }

  FILE *output_file = fopen("../output/output.dat", "w");
  if (!output_file) {
    printf("Errore inizializzazione file\n");
    return 1;
  }

  if (1) {
    printf("ErfcTable generation started\n");
    FILE *tableErfcFile = fopen("../output/erfc_table.dat", "w");
    if (!tableErfcFile) {
      printf("Errore inizializzazione file\n");
      return 1;
    }
    generateErfcTable(ERFC_TABLE_IN, ERFC_TABLE_FIN, ERFC_TABLE_PRECISION, tableErfcFile);
    printf("ErfcTable generation completetd\n");
  }

  // NOTE - erfc speed test

  clock_t start = clock();
  for (size_t i = 0; i < N_PARTICLES; i++) {
    for (size_t j = 0; j < N_PARTICLES; j++) {
      double x = exp(i + j);
    }
  }
  clock_t end = clock();
  double time_spent = (double)(end - start) / CLOCKS_PER_SEC;
  printf("Erf time: %lf ms\n", time_spent * 1000);

  const double TIME_IN = 0;
  const double TIME_END = 1;
  const double DELTA_T = 1e-3;
  const int N_STEPS = (TIME_END - TIME_IN) / DELTA_T;

  // BRUTE FORCE INIT ------------
  // NOTE: non tiene conto della possibilità di due particelle sovrapposte
  for (size_t i = 0; i < N_PARTICLES * 3; i += 3) {
    pos_array[i + 0] = rand() / (RAND_MAX + 1.0) * CELL_L;
    pos_array[i + 1] = rand() / (RAND_MAX + 1.0) * CELL_L;
    pos_array[i + 2] = rand() / (RAND_MAX + 1.0) * CELL_L;
  }

  // -----------------------------

  for (size_t i = 0; i < N_STEPS; i++) {
    double cur_t = TIME_IN + i * DELTA_T;
    clock_t start = clock();
    int result = verletPropagationStep(pos_array, vel_array, forces_array_ptr, masses_array, N_PARTICLES, cur_t, DELTA_T, edwald_summation, NULL, 0);
    if (result) {
      printf("Errore verlet propagation:\nIndice: %d\nCur_time: %lf\n", i, cur_t);
      return 1;
    }

    clock_t end = clock();
    double time_spent = (double)(end - start);
    printf("Verlet time: %lf ms\n", time_spent);

    // printf("ENERGIA TOTALE: %lf a tempo %lf\n", kinetic_energy(vel_array, masses_array, N_PARTICLES) + coulomb_potential_energy(pos_array, N_PARTICLES), cur_t);
  }

  free(pos_array);
  free(vel_array);
  free(masses_array);
  free(*forces_array_ptr);
  free(forces_array_ptr);
  return 0;
}
