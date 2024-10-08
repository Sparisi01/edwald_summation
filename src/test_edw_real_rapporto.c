#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <tgmath.h>

#include "includes/coulomb/coulomb_potential.h"
#include "includes/ewald/edwald.h"
#include "includes/utils/statistic.h"

int _N_PARTICLES = 64;
double _DENSITY = 0.01;
double _CELL_LENGHT = 1;
double _SIGMA_VELOCITIES = 1.;

int writeParticlesPositions(Particle *particles, int n_particles, FILE *file)
{
    if (!file) return 1;
    if (!particles) return 1;

    for (size_t i = 0; i < n_particles; i++)
    {
        fprintf(file, "%.5E;%.5E;%.5E\n", particles[i].x, particles[i].y, particles[i].z);
    }
    return 0;
}

int main(int argc, char const *argv[])
{

    printf("%lf\n", alpha_by_precision(1e-3));

    printlnElementSymbol(107);
    srand(RAND_SEED);
    System system;
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
        system.particles[i].charge = (i % 2 == 0) ? 1 : -1;
        charge_sum += system.particles[i].charge;
    }

    printf("Total charge: %lf\n", charge_sum);

    _K_RANGE_EWALD = 25;

    double ALPHA_MAX = 15;
    double ALPHA_MIN = 2;
    _CUTOFF = 5*(_CELL_LENGHT/2);
    int N_ALPHA = 29;

    FILE *file_comparison_real_rec = fopen("../src/data/comparison_real_rec.csv", "w");

    for (size_t i = 0; i < N_ALPHA; i++)
    {
        double tmp_alpha = ALPHA_MIN + (ALPHA_MAX - ALPHA_MIN) / ((double)N_ALPHA) * i;
        double real = real_space_coulomb_energy(&system, tmp_alpha) / _N_PARTICLES;
        double rec = reciprocal_space_coulomb_energy(&system, tmp_alpha) / _N_PARTICLES;
        double self = self_coulomb_energy(&system, tmp_alpha) / _N_PARTICLES;

        printf("%.5E;%.5E;%.5E;%.5E\n", tmp_alpha, real, rec, self);
        fprintf(file_comparison_real_rec, "%.16E;%.16E;%.16E;%.16E;%.16E\n", tmp_alpha, real, rec, self, real + rec - self);
    }
}
