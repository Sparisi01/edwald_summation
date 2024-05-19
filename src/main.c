#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <tgmath.h>

#include "includes/complex_utils.h"
#include "includes/coulomb/coulomb_potential.h"
#include "includes/ewald/edwald.h"
#include "includes/statistic.h"

int _N_PARTICLES = 100;
double _DENSITY = 2;
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
    srand(RAND_SEED);
    System system;
    system.n_particles = _N_PARTICLES / _DENSITY;
    system.cell_lenght = _CELL_LENGHT;
    system.particles = (Particle *)malloc(sizeof(Particle) * system.n_particles);
    if (!system.particles)
    {
        perror("Error, malloc 'system.particles' returned NULL:");
        exit(EXIT_FAILURE);
    }
    // INITIALIZATION ------------
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
    FILE *file2 = fopen("./data/particel_start_pos.csv", "w");
    writeParticlesPositions(system.particles, system.n_particles, file2);
    // DO THINGS
    FILE *file = fopen("./data/range_variabile_3.csv", "w");
    _CUTOFF = _CELL_LENGHT / 2;
    _ALPHA = 3 / _CUTOFF;

    // Stima errore se particelle disposte casualmente

    double max_error_short = erfc(_ALPHA * _CUTOFF) / _CUTOFF;
    double volume_esterno_sfera = pow(_CELL_LENGHT, 3) - 4. / 3. * PI * pow(_CELL_LENGHT / 2., 3);
    double n_interazioni = system.n_particles * (system.n_particles - 1) / 2.;
    printf("Stima errore totale: %.5E\n", max_error_short * n_interazioni * volume_esterno_sfera);
    printf("Stima errore singolo: %.5E\n", max_error_short);

    for (size_t i = 0; i < 20; i++)
    {
        _K_RANGE = i;
        double pot = getEdwaldPotential(&system);
        // fprintf(file, "%d;%lf\n", i, pot);

        printf("ED: %.5E\n", pot);
    }

    for (size_t i = 0; i < 10; i++)
    {
        _R_RANGE = i;
        // double pot = getCoulombPotential(&system);
        double pot = getCoulombPotential(&system);
        fprintf(file, "%d;%lf\n", i, pot);
        printf("C: %.5E\n", pot);
    }
}
