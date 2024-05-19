#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <tgmath.h>

#include "includes/complex_utils.h"
#include "includes/ewald/edwald_method.h"
#include "includes/statistic.h"
#include <math.h>

int _N_PARTICLES = 1000;
double _DENSITY = 1;
double _CELL_LENGHT = 2;
double _SIGMA_VELOCITIES = 1.;

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

    // DO THINGS
    FILE *file = fopen("./data/range_variabile_3.csv", "w");
    //_CUTOFF = _CELL_LENGHT / 2;
    _ALPHA = 5;

    // Stima errore se particelle disposte casualmente

    double max_error_short = erfc(_ALPHA * _CUTOFF) / _CUTOFF;
    double volume_esterno_sfera = pow(_CELL_LENGHT, 3) - 4. / 3. * PI * pow(_CELL_LENGHT / 2., 3);
    double n_interazioni = system.n_particles * (system.n_particles - 1) / 2.;
    printf("Stima errore short: %.5E\n", max_error_short * n_interazioni * volume_esterno_sfera);

    for (size_t i = 0; i < 25; i++)
    {
        double pot = getEdwaldPotential(&system);
        fprintf(file, "%d;%lf\n", i, pot);
        _K_RANGE = i;
        printf("%.5E\n", pot);
    }
}
