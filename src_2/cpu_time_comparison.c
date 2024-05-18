#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "constants.h"
#include "coulomb.h"
#include "edwald_summation.h"
#include "statistic.h"

int main(int argc, char const *argv[])
{
    FILE *file = fopen("../output/cpu_time_edwald_parallelized_10.dat", "w");

    int N[11] = {50, 70, 100, 150, 200, 300, 400, 500, 600, 700, 1000};

    for (size_t i = 0; i < 11; i++)
    {

        System system;
        system.n_particles = N[i];
        system.particles = (Particle *)malloc(sizeof(Particle) * system.n_particles);
        if (!system.particles)
        {
            perror("Error, malloc 'system.particles' returned NULL:");
            exit(EXIT_FAILURE);
        }

        for (size_t i = 0; i < system.n_particles; i++)
        {
            system.particles[i].x = randUnif(-CELL_LENGHT / 2, CELL_LENGHT / 2);
            system.particles[i].y = randUnif(-CELL_LENGHT / 2, CELL_LENGHT / 2);
            system.particles[i].z = randUnif(-CELL_LENGHT / 2, CELL_LENGHT / 2);

            system.particles[i].vx = randGauss(0, SIGMA_VELOCITIES);
            system.particles[i].vy = randGauss(0, SIGMA_VELOCITIES);
            system.particles[i].vz = randGauss(0, SIGMA_VELOCITIES);

            system.particles[i].mass = 1;
            system.particles[i].charge = 0.1;

            // printf("%1.10E\n", system.particles[i].vz);
        }

        system.forces = edwald_summation(&system, NULL);

        long sum = 0;
        int N_steps = 1000;
        for (size_t i = 0; i < N_steps; i++)
        {
            clock_t time_start = clock();

            Vec3 *result = edwald_summation_parallelized(&system, NULL);

            clock_t time_end = clock();

            sum += time_end - time_start;
        }

        double cpu_time = (double)sum / CLOCKS_PER_SEC / N_steps;
        printf("CPU time for cicle: %lf s\n", cpu_time);
        printf("--------------------\n");

        fprintf(file, "%d %lf\n", N[i], cpu_time);
    }

    return 0;
}
