#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
// commento
#include "constants.h"
#include "edwald_summation.h"
#include "observables.h"
#include "statistic.h"
#include "structures.h"
#include "verlet_propagation.h"

int writeParticlesPositions(Particle *particles, int n_particles, FILE *file)
{
    if (!file) return 1;
    if (!particles) return 1;

    for (size_t i = 0; i < n_particles; i++)
    {
        fprintf(file, "%.5E %.5E %.5E\n", particles[i].x, particles[i].y, particles[i].z);
    }
    return 0;
}

int main(int argc, char const *argv[])
{
    srand(RAND_SEED);

    FILE *file_start_particles_pos = fopen("../output/start_pos.dat", "w");
    if (!file_start_particles_pos)
    {
        perror("Error trying to open 'file_start_particles_pos':");
        exit(EXIT_FAILURE);
    }

    FILE *file_end_particles_pos = fopen("../output/end_pos.dat", "w");
    if (!file_end_particles_pos)
    {
        perror("Error trying to open 'file_end_particles_pos':");
        exit(EXIT_FAILURE);
    }

    FILE *file_thermodinamic = fopen("../output/energy.dat", "w");
    if (!file_thermodinamic)
    {
        perror("Error trying to open 'file_thermodinamic':");
        exit(EXIT_FAILURE);
    }

    struct System system;
    system.n_particles = N_PARTICLES;
    system.particles = (struct Particle *)malloc(sizeof(struct Particle) * system.n_particles);
    if (!system.particles)
    {
        perror("Error, malloc 'system.particles' returned NULL:");
        exit(EXIT_FAILURE);
    }

    //----------------------------------------------
    double start_time = 0;
    double end_time = 6;
    double time_step = 1e-3;

    int n_time_step = (end_time - start_time) / time_step;

    struct Observables observables;

    observables.kinetic_energy = (double *)calloc(n_time_step, sizeof(double));
    if (!observables.kinetic_energy)
    {
        perror("Error, calloc 'observables.kinetic_energy' returned NULL:");
        exit(EXIT_FAILURE);
    }

    observables.potential_energy = (double *)calloc(n_time_step, sizeof(double));
    if (!observables.potential_energy)
    {
        perror("Error, calloc 'observables.potential_energy' returned NULL:");
        exit(EXIT_FAILURE);
    }

    observables.pressure = (double *)calloc(n_time_step, sizeof(double));
    if (!observables.pressure)
    {
        perror("Error, calloc 'observables.pressure' returned NULL:");
        exit(EXIT_FAILURE);
    }

    observables.temperature = (double *)calloc(n_time_step, sizeof(double));
    if (!observables.temperature)
    {
        perror("Error, calloc 'observables.temperature' returned NULL:");
        exit(EXIT_FAILURE);
    }

    observables.time = (double *)calloc(n_time_step, sizeof(double));
    if (!observables.time)
    {
        perror("Error, calloc 'observables.time' returned NULL:");
        exit(EXIT_FAILURE);
    }

    // INITIALIZATION ------------

    for (size_t i = 0; i < system.n_particles; i++)
    {
        system.particles[i].x = rand() / (RAND_MAX + 1.0) * CELL_LENGHT - CELL_LENGHT / 2;
        system.particles[i].y = rand() / (RAND_MAX + 1.0) * CELL_LENGHT - CELL_LENGHT / 2;
        system.particles[i].z = rand() / (RAND_MAX + 1.0) * CELL_LENGHT - CELL_LENGHT / 2;

        system.particles[i].vx = randGauss(0, SIGMA_VELOCITIES);
        system.particles[i].vy = randGauss(0, SIGMA_VELOCITIES);
        system.particles[i].vz = randGauss(0, SIGMA_VELOCITIES);

        system.particles[i].mass = 1;
        system.particles[i].charge = 0.1;
    }

    printf("--------------------\n");
    printf("AVVIO SIMULAZIONE\n");
    printf("--------------------\n");
    printf("N Particelle: %d\n", system.n_particles);
    printf("Volume: %lf\n", CELL_VOLUME);
    printf("Δt, dt: %lf, %lf\n", end_time - start_time, time_step);
    printf("--------------------\n");

    system.forces = edwald_summation(&system, NULL);

    writeParticlesPositions(system.particles, system.n_particles, file_start_particles_pos);

    int n_time_measure = 10;
    clock_t comput_time_start = 0;
    clock_t sum_comput_times = 0;

    for (size_t i = 0; i < n_time_step; i++)
    {

        comput_time_start = clock();

        system.time = start_time + i * time_step;

        int result = verletPropagationStep(&system, time_step, edwald_summation, NULL);
        if (result)
        {
            perror("Error in verlet propagation step: ");
            exit(EXIT_FAILURE);
        }

        observables.time[i] = system.time;
        observables.kinetic_energy[i] = kinetic_energy(system.particles, system.n_particles);
        observables.potential_energy[i] = potential_energy(system.particles, system.n_particles);
        observables.temperature[i] = 2. / 3. * observables.kinetic_energy[i] / system.n_particles;
        observables.pressure[i] = 0;

        // Exclude the first one because it contains the edwald table generation
        if (i < n_time_measure && i != 0)
        {
            sum_comput_times += clock() - comput_time_start;
        }
        if (i == n_time_measure)
        {
            double time_spent = (double)sum_comput_times / n_time_measure;
            int tot_time_sec = time_spent * n_time_step / CLOCKS_PER_SEC;
            int sec = tot_time_sec % 60;        // Estimated time seconds
            int min = (tot_time_sec / 60) % 60; // Estimated time minutes
            int h = tot_time_sec / 3600;        // Estimated time hours
            printf("--------------------\nVerlet time: %d ms\nTotale estimated time: %d:%d:%d [h:m:s]\n--------------------\n", (int)time_spent, h, min, sec);
        }
    }

    writeParticlesPositions(system.particles, system.n_particles, file_end_particles_pos);

    // Save energies in file and print statistics
    for (size_t i = 0; i < n_time_step; i++)
    {
        fprintf(file_thermodinamic, "%.10E %.10E %.10E %.10E %.10E %.10E\n",
                start_time + i * time_step,
                observables.kinetic_energy[i] + observables.potential_energy[i],
                observables.kinetic_energy[i],
                observables.potential_energy[i], observables.temperature[i],
                observables.pressure[i]);
    }

    int equilibrium_step = 0;

    double mean_temperature = mean(observables.temperature, equilibrium_step, n_time_step);
    double std_temperature = stddev(observables.temperature, equilibrium_step, n_time_step);
    printf("TEMPERATURE = %.5E ± %.5E\n", mean_temperature, std_temperature);

    fclose(file_start_particles_pos);
    fclose(file_end_particles_pos);

    free(system.particles);
    free(observables.kinetic_energy);
    free(observables.potential_energy);
    free(observables.pressure);
    free(observables.temperature);
    free(observables.time);

    exit(EXIT_SUCCESS);
}
