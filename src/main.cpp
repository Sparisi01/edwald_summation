#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "edwald_summation.h"
#include "verlet_propagation.h"

using namespace std;

#define N_PARTICLE 1

double* forza_elastica(double t, double* pos, double* vel, int n_particles, double* args, int n_args) {
    double* new_forces = (double*)calloc(n_particles * 3, sizeof(double));
    if (!new_forces) {
        return NULL;
    }

    for (size_t i = 0; i < n_particles; i += 3) {
        new_forces[i + 0] = -5 * pos[i];
        new_forces[i + 1] = 0;
        new_forces[i + 2] = 0;
    }

    return new_forces;
}

int main(int argc, char const* argv[]) {
    double* pos_array = (double*)calloc(N_PARTICLE * 3, sizeof(double));
    double* vel_array = (double*)calloc(N_PARTICLE * 3, sizeof(double));
    double** forces_array_ptr = (double**)malloc(sizeof(double*));
    *forces_array_ptr = (double*)calloc(N_PARTICLE * 3, sizeof(double));
    double* masses_array = (double*)calloc(N_PARTICLE, sizeof(double));

    if (!pos_array || !vel_array || !forces_array_ptr || !*forces_array_ptr || !masses_array) {
        printf("Errore inizializzazione array\n");
        return 1;
    }

    for (size_t i = 0; i < N_PARTICLE; i++) {
        masses_array[i] = 1;
    }

    FILE* output_file = fopen("../output/output.dat", "w");
    if (!output_file) {
        printf("Errore inizializzazione file\n");
        return 1;
    }

    const double TIME_IN = 0;
    const double TIME_END = 10;
    const double DELTA_T = 1e-3;
    const int N_STEPS = (TIME_END - TIME_IN) / DELTA_T;

    // Init
    pos_array[0] = 50;
    for (size_t i = 0; i < N_STEPS; i++) {
        double cur_t = TIME_IN + i * DELTA_T;
        int result = verletPropagationStep(pos_array, vel_array, forces_array_ptr, masses_array, N_PARTICLE, cur_t, DELTA_T, forza_elastica);
        if (result) {
            printf("Errore verlet propagation:\nIndice: %d\nCur_time: %lf\n", i, cur_t);
            return 1;
        }
        fprintf(output_file, "%lf %lf %lf %lf\n", cur_t, pos_array[0], pos_array[1], pos_array[2]);
    }

    free(pos_array);
    free(vel_array);
    free(masses_array);
    free(*forces_array_ptr);
    free(forces_array_ptr);
    return 0;
}
