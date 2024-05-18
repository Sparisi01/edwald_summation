#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "constants.h"
#include "coulomb.h"
#include "edwald_summation.h"
#include "statistic.h"

int main(int argc, char const *argv[])
{
    FILE *file = fopen("../output/convergenze/convergenza_coulomb.csv", "w");
    FILE *file2 = fopen("../output/convergenze/convergenza_ewald.csv", "w");

    System system;
    system.n_particles = N_PARTICLES;
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

        system.particles[i].charge = (i % 2 == 0) ? 0.1 : -0.1; // Total charge null if N even
    }

    double charge_sum = 0;
    for (size_t i = 0; i < system.n_particles; i++)
    {
        charge_sum += system.particles[i].charge;
    }
    
    printf("TOTAL CHARGE: %lf\n", charge_sum);

    /* fprintf(file, "RANGE;POTENTIAL\n");
    for (size_t i = 0; i < 20; i++)
    {
        COULOMB_CELL_RANGE = i;
        double coulomb = coulomb_potential(&system);
        fprintf(file, "%d;%.8E\n", i, coulomb);
        printf("%d\n", i);
    } */

    double ewald = edwald_potential(&system);
    printf("%.8E\n", ewald);
    return 0;
}
