#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <tgmath.h>

#include "includes/complex_utils.h"
#include "includes/ewald/edwald_method.h"
#include "includes/statistic.h"

double _CELL_LENGHT = 2.;
double _SIGMA_VELOCITIES = 1.;

int main(int argc, char const *argv[])
{
    System system;
    system.n_particles = 10;
    system.cell_lenght = _CELL_LENGHT;
    system.particles = (Particle *)malloc(sizeof(Particle) * system.n_particles);
    if (!system.particles)
    {
        perror("Error, malloc 'system.particles' returned NULL:");
        exit(EXIT_FAILURE);
    }
    // INITIALIZATION ------------

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
    }

    setAlphaByError(1E-3);
    double pot = getEdwaldPotential(&system);
    printf("%.5E", pot);
}
