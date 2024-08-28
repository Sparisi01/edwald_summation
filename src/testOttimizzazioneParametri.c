#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <tgmath.h>
#include <time.h>

#include "includes/coulomb/coulomb_potential.h"
#include "includes/ewald/edwald.h"
#include "includes/utils/statistic.h"

int _N_PARTICLES = 100;
double _CELL_LENGHT = 1;
double _SIGMA_VELOCITIES = 1.;
double _CHARGE = 1.;

int main(int argc, char const *argv[])
{
    _ALPHA = 9;
    _CUTOFF = 0.5;
    _K_RANGE_EWALD = 12;

    FILE *tempi_file = fopen("data/execution_times_file.csv", "w");

    System system;
    double pot2 = 0;

    for (size_t j = 5; j < 100; j++)
    {

        _N_PARTICLES = 100 * j;

        system.n_particles = _N_PARTICLES;
        system.cell_lenght = _CELL_LENGHT;
        system.particles = (Particle *)malloc(sizeof(Particle) * system.n_particles);
        if (!system.particles)
        {
            perror("Error, malloc 'system.particles' returned NULL:");
            exit(EXIT_FAILURE);
        }

        //==>INITIALIZATION<==//
        double charge_sum = 0;
        for (size_t i = 0; i < system.n_particles; i++)
        {
            system.particles[i].x = randUnif(-system.cell_lenght / 2, system.cell_lenght / 2);
            system.particles[i].y = randUnif(-system.cell_lenght / 2, system.cell_lenght / 2);
            system.particles[i].z = randUnif(-system.cell_lenght / 2, system.cell_lenght / 2);

            system.particles[i].vx = randGauss(0, _SIGMA_VELOCITIES);
            system.particles[i].vy = randGauss(0, _SIGMA_VELOCITIES);
            system.particles[i].vz = randGauss(0, _SIGMA_VELOCITIES);

            system.particles[i].mass = 1;
            system.particles[i].charge = (i % 2 == 0) ? 0.1 : -0.1;
            charge_sum += system.particles[i].charge;
        }

        printf("Total charge: %lf\n", charge_sum);

        double Q = 1. * system.n_particles;

        optimizeParameter(1e-3, _CELL_LENGHT, system.n_particles, Q);

        clock_t tic = clock();

        double short_range = real_space_coulomb_energy(&system, _ALPHA);
        double long_range = reciprocal_space_coulomb_energy(&system, _ALPHA);
        double self = self_coulomb_energy(&system, _ALPHA);

        clock_t toc = clock();

        double pot = (short_range + long_range - self);

        // PLOT CORRETTO
        /* _ALPHA = 9;
        _CUTOFF = 2 * (_CELL_LENGHT / 2);
        _K_RANGE_EWALD = 20;

        short_range = real_space_coulomb_energy(&system, _ALPHA);
        long_range = reciprocal_space_coulomb_energy(&system, _ALPHA);
        self = self_coulomb_energy(&system, _ALPHA);

        double pot2 = (short_range + long_range - self); */

        printf("{N: %d, Alpha: %.3E, nr: %d, nk: %d, Real:%.3E, Rec:%.3E, Self:%.3E, TotEnergy:%.10E}\n", system.n_particles, _ALPHA, _R_RANGE_EWALD, _K_RANGE_EWALD, short_range, long_range, self, pot);
        fprintf(tempi_file, "%d;%f;%d;%f;%.5E\n", system.n_particles, _CUTOFF, _K_RANGE_EWALD, (double)(toc - tic), fabs(pot - pot2));
    }

    return 0;
}
