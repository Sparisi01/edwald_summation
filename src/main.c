#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "constants.h"
#include "edwald_summation.h"
#include "statistic.h"
#include "thermodinamics.h"
#include "verlet_propagation.h"

/**
 * TODO: main todo
 * [] ricerca sigma sensato.
 * [x] verifica sistemi noti.
 * [x] Cutof potenziale.
 * [] Verificare se è necessario anche per la somma nello spazio
 * reciproco considerare la coppia di particelle più vicina.
 */

int generateErfcTable(double in, double fin, double precision, FILE *tableFile)
{
    if (!tableFile)
    {
        return 1;
    }

    int N = (fin - in) / precision;
    for (size_t i = 0; i < N; i++)
    {
        fprintf(tableFile, "%lf\n", erfc(in + i * precision));
    }

    return 0;
}

double *loadErfcTable(double in, double fin, double precision, FILE *tableFile)
{
    int N = (fin - in) / precision;
    double *array = (double *)malloc(sizeof(double) * N);
    if (!array)
    {
        return NULL;
    }
    if (!tableFile)
    {
        return NULL;
    }
    else
    {
        int m = 0;
        double next;
        do
        {
            next = fscanf(tableFile, "%lf", &array[m++]);
        } while (m != N);
    }
    return array;
}

int saveParticelsPositions(double *pos, int n_particles, FILE *file)
{
    if (!file) return 1;

    for (size_t i = 0; i < n_particles * 3; i += 3)
    {
        fprintf(file, "%lf %lf %lf\n", pos[i + 0], pos[i + 1], pos[i + 2]);
    }
    return 0;
}

int main(int argc, char const *argv[])
{
    // Set Seed for rnd numbers
    srand(SEED);

    double *pos_array = (double *)calloc(N_PARTICLES * 3, sizeof(double));
    double *vel_array = (double *)calloc(N_PARTICLES * 3, sizeof(double));
    double **forces_array_ptr = (double **)malloc(sizeof(double *));
    *forces_array_ptr = (double *)calloc(N_PARTICLES * 3, sizeof(double));
    double *masses_array = (double *)calloc(N_PARTICLES, sizeof(double));
    double *force_charge_array = (double *)calloc(N_PARTICLES, sizeof(double));

    if (!pos_array || !vel_array || !forces_array_ptr || !*forces_array_ptr || !masses_array || !force_charge_array)
    {
        perror("Error init particles properties arrays: ");
        exit(EXIT_FAILURE);
    }

    FILE *start_pos_file = fopen("../output/start_pos.dat", "w");
    FILE *end_pos_file = fopen("../output/end_pos.dat", "w");
    if (!start_pos_file || !end_pos_file)
    {
        perror("Error start-end pos file: ");
        exit(EXIT_FAILURE);
    }

    /* if (0) {
        printf("ErfcTable generation started\n");
        FILE* tableErfcFile = fopen("../output/erfc_table.dat", "w");
        if (!tableErfcFile) {
            printf("Errore inizializzazione file\n");
            exit(EXIT_FAILURE);
        }
        generateErfcTable(ERFC_TABLE_IN, ERFC_TABLE_FIN, ERFC_TABLE_PRECISION, tableErfcFile);
        printf("ErfcTable generation completetd\n");
        fclose(tableErfcFile);
    } */

    // Load ERFC TABLE
    /* FILE* tableErfcFile = fopen("../output/erfc_table.dat", "r");
    if (!tableErfcFile) exit(EXIT_FAILURE);

    double* erfcTable = loadErfcTable(ERFC_TABLE_IN, ERFC_TABLE_FIN, ERFC_TABLE_PRECISION, tableErfcFile);
    if (!erfcTable) {
        printf("Error erfc table\n");
        exit(EXIT_FAILURE);
    }
    fclose(tableErfcFile);
    printf("Erfc table loaded\n"); */

    // Energy arrays
    FILE *file_energy = fopen("../output/energy.dat", "w");
    double *kinetic_energy_array = (double *)malloc(N_STEPS * sizeof(double));
    double *potential_energy_array = (double *)malloc(N_STEPS * sizeof(double));
    double *energy_array = (double *)malloc(N_STEPS * sizeof(double));
    if (!potential_energy_array || !kinetic_energy_array || !energy_array || !file_energy)
    {
        perror("Error energy array init: ");
        exit(EXIT_FAILURE);
    }

    // INITIALIZATION ------------
    // NOTE: non tiene conto della possibilità di due particelle sovrapposte
    for (size_t i = 0; i < N_PARTICLES * 3; i += 3)
    {
        pos_array[i + 0] = rand() / (RAND_MAX + 1.0) * CELL_L - CELL_L / 2;
        pos_array[i + 1] = rand() / (RAND_MAX + 1.0) * CELL_L - CELL_L / 2;
        pos_array[i + 2] = rand() / (RAND_MAX + 1.0) * CELL_L - CELL_L / 2;

        vel_array[i + 0] = 0;
        vel_array[i + 1] = 0;
        vel_array[i + 2] = 0;

        force_charge_array[i / 3] = Q;

        masses_array[i / 3] = 1;
    }

    /* pos_array[0 + 0] = 0.3;
    pos_array[3 + 0] = -0.3;
    pos_array[6 + 1] = 0.3;
    pos_array[9 + 1] = -0.3;
    force_charge_array[0] = Q;
    force_charge_array[1] = Q;
    force_charge_array[2] = Q;
    force_charge_array[3] = Q; */

    // -----------------------------

    printf("--------------------\n");
    printf("AVVIO SIMULAZIONE\n");
    printf("--------------------\n");
    printf("N Particelle: %d\n", (int)N_PARTICLES);
    printf("Dimensione cella: %lf\n", (double)CELL_L);
    printf("Δt, dt: %lf, %lf\n", (double)(TIME_END - TIME_IN), (double)DELTA_T);
    printf("--------------------\n");

    saveParticelsPositions(pos_array, N_PARTICLES, start_pos_file); // Save particles starting position

    const int count_dow_time_measure = 5;
    clock_t start = 0;
    clock_t end = 0;
    for (size_t i = 0; i < N_STEPS; i++)
    {
        double cur_t = TIME_IN + i * DELTA_T;

        /** At the first count_dow_time_measure step measure time of a verletPropagationStep.
         * Multiply this time for N_STEPS in order to estimate the total
         * computational time.
         */
        if (i < count_dow_time_measure) start += clock();

        int result = verletPropagationStep(pos_array, vel_array, forces_array_ptr, masses_array, N_PARTICLES, cur_t, DELTA_T, edwald_summation, (void *)force_charge_array);
        if (result)
        {
            perror("Error in verlet propagation step: ");
            exit(EXIT_FAILURE);
        }

        kinetic_energy_array[i] = kinetic_energy(vel_array, masses_array, N_PARTICLES);
        potential_energy_array[i] = potential_energy(pos_array, force_charge_array, N_PARTICLES);
        energy_array[i] = kinetic_energy_array[i] + potential_energy_array[i];

        if (i < count_dow_time_measure) end += clock();
        if (i == count_dow_time_measure - 1)
        {
            double time_spent = (double)(end - start) / count_dow_time_measure;
            int tot_time_sec = time_spent * N_STEPS / CLOCKS_PER_SEC;
            int sec = tot_time_sec % 60;        // Estimated time seconds
            int min = (tot_time_sec / 60) % 60; // Estimated time minutes
            int h = tot_time_sec / 3600;        // Estimated time hours
            printf("--------------------\nVerlet time: %d ms\nTempo totale stimato: %d:%d:%d [h:m:s]\n--------------------\n", (int)time_spent, h, min, sec);
        }
    }

    saveParticelsPositions(pos_array, N_PARTICLES, end_pos_file); // Save particles final position

    // Save energies in file and print statistics
    for (size_t i = 0; i < N_STEPS; i++)
    {
        fprintf(file_energy, "%.10E %.10E %.10E %.10E\n", DELTA_T * i + TIME_IN, energy_array[i], potential_energy_array[i], kinetic_energy_array[i]);
    }

    printf("ENERGY: %.5E ± %.5E \n", mean(energy_array, 0, N_STEPS), stddev(energy_array, 0, N_STEPS));

    /**
     * - Free all memory
     * - Close all file stream
     */
    fclose(start_pos_file);
    fclose(end_pos_file);
    // free(erfcTable);
    free(pos_array);
    free(vel_array);
    free(masses_array);
    free(*forces_array_ptr);
    free(forces_array_ptr);
    free(force_charge_array);
    free(kinetic_energy_array);
    free(potential_energy_array);
    exit(EXIT_SUCCESS);
}
