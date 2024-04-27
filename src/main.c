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
 * [] RICERCA SIGMA SENSATO.
 * [] VERIFICA SISTEMI NOTI.
 * [x] Cutof potenziale.
 * [] Verificare se è necessario anche per la somma nello spazio
 * reciproco considerare la coppia di particelle più vicina.
 */

int generateErfcTable(double in, double fin, double precision, FILE* tableFile) {
    if (!tableFile) {
        return 1;
    }

    int N = (fin - in) / precision;
    for (size_t i = 0; i < N; i++) {
        fprintf(tableFile, "%lf\n", erfc(in + i * precision));
    }

    return 0;
}

double* loadErfcTable(double in, double fin, double precision, FILE* tableFile) {
    int N = (fin - in) / precision;
    double* array = (double*)malloc(sizeof(double) * N);
    if (!array) {
        return NULL;
    }
    if (!tableFile) {
        return NULL;
    } else {
        int m = 0;
        double next;
        do {
            next = fscanf(tableFile, "%lf", &array[m++]);
        } while (m != N);
    }
    return array;
}

int main(int argc, char const* argv[]) {
    // Set Seed for rnd numbers
    srand(SEED);

    double* pos_array = (double*)calloc(N_PARTICLES * 3, sizeof(double));
    double* vel_array = (double*)calloc(N_PARTICLES * 3, sizeof(double));
    double** forces_array_ptr = (double**)malloc(sizeof(double*));
    *forces_array_ptr = (double*)calloc(N_PARTICLES * 3, sizeof(double));
    double* masses_array = (double*)calloc(N_PARTICLES, sizeof(double));
    double* force_charge_array = (double*)calloc(N_PARTICLES, sizeof(double));

    if (!pos_array || !vel_array || !forces_array_ptr || !*forces_array_ptr || !masses_array || !force_charge_array) {
        printf("Errore inizializzazione array\n");
        exit(EXIT_FAILURE);
    }

    for (size_t i = 0; i < N_PARTICLES; i++) {
        masses_array[i] = 1;
    }

    FILE* output_file = fopen("../output/output.dat", "w");
    if (!output_file) {
        printf("Errore inizializzazione file\n");
        exit(EXIT_FAILURE);
    }

    if (0) {
        printf("ErfcTable generation started\n");
        FILE* tableErfcFile = fopen("../output/erfc_table.dat", "w");
        if (!tableErfcFile) {
            printf("Errore inizializzazione file\n");
            exit(EXIT_FAILURE);
        }
        generateErfcTable(ERFC_TABLE_IN, ERFC_TABLE_FIN, ERFC_TABLE_PRECISION, tableErfcFile);
        printf("ErfcTable generation completetd\n");
        fclose(tableErfcFile);
    }

    // Load ERFC TABLE
    FILE* tableErfcFile = fopen("../output/erfc_table.dat", "r");

    double* erfcTable = loadErfcTable(ERFC_TABLE_IN, ERFC_TABLE_FIN, ERFC_TABLE_PRECISION, tableErfcFile);
    if (!erfcTable) {
        exit(EXIT_FAILURE);
    }
    fclose(tableErfcFile);
    printf("Erfc table loaded");

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
        force_charge_array[i / 3] = Q;
    }

    // -----------------------------

    for (size_t i = 0; i < N_STEPS; i++) {
        double cur_t = TIME_IN + i * DELTA_T;
        clock_t start = clock();
        int result = verletPropagationStep(pos_array, vel_array, forces_array_ptr, masses_array, N_PARTICLES, cur_t, DELTA_T, edwald_summation, (void*)force_charge_array);
        // int result = verletPropagationStep(pos_array, vel_array, forces_array_ptr, masses_array, N_PARTICLES, cur_t, DELTA_T, edwald_summation_table, erfcTable);
        if (result) {
            printf("Errore verlet propagation:\nIndice: %d\nCur_time: %lf\n", i, cur_t);
            exit(EXIT_FAILURE);
        }

        clock_t end = clock();
        double time_spent = (double)(end - start);
        printf("Verlet time: %lf ms\n", time_spent);

        printf("ENERGIA TOTALE: %lf a tempo %lf\n", kinetic_energy(vel_array, masses_array, N_PARTICLES) + coulomb_potential_energy(pos_array, N_PARTICLES), cur_t);
        printf("Posizione: %lf %lf %lf \n", pos_array[0], pos_array[1], pos_array[2]);
    }

    free(erfcTable);
    free(pos_array);
    free(vel_array);
    free(masses_array);
    free(*forces_array_ptr);
    free(forces_array_ptr);
    free(force_charge_array);
    exit(EXIT_SUCCESS);
}
