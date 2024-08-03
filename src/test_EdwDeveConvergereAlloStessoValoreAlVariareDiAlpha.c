#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <tgmath.h>

#include "includes/coulomb/coulomb_potential.h"
#include "includes/ewald/edwald.h"
#include "includes/utils/statistic.h"

int _N_PARTICLES = 10;
double _DENSITY = 0.01;
double _CELL_LENGHT = 1;
double _SIGMA_VELOCITIES = 1.;
int _MAX_RANGE = 50;

typedef struct alpha_file
{
    double alpha;
    char *file;
} alpha_file;

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
    FILE *file_start_position = fopen("./data/particel_start_pos.csv", "w");

    writeParticlesPositions(system.particles, system.n_particles, file_start_position);

    double last_pot;

    _CUTOFF = _CELL_LENGHT / 2;

    alpha_file config[7] = {
        {
            .alpha = 1.,
            .file = "./data/convergenza_range/range_variabile_edw_1.csv",

        },
        {
            .alpha = 2.,
            .file = "./data/convergenza_range/range_variabile_edw_2.csv",

        },
        {
            .alpha = 3.5,
            .file = "./data/convergenza_range/range_variabile_edw_3_5.csv",

        },
        {
            .alpha = 5.,
            .file = "./data/convergenza_range/range_variabile_edw_5.csv",

        },
        {
            .alpha = 7.,
            .file = "./data/convergenza_range/range_variabile_edw_7.csv",

        },
        {
            .alpha = 9.,
            .file = "./data/convergenza_range/range_variabile_edw_9.csv",
        },
        {
            .alpha = 15.,
            .file = "./data/convergenza_range/range_variabile_edw_15.csv",
        },
    };

    for (size_t i = 0; i < 7; i++)
    {
        FILE *file_convergenza_edw = fopen(config[i].file, "w");
        if (!file_convergenza_edw) exit(EXIT_FAILURE);
        _ALPHA = config[i].alpha / _CUTOFF;
        last_pot = 0;
        printf("---------------------\n");
        printf("ALPHA: %lf\n", _ALPHA);
        for (size_t j = 0; j < _MAX_RANGE; j++)
        {
            _K_RANGE = j;
            double pot = ewald_energy(&system);
            double rel_error = fabs((last_pot - pot) / pot);

            fprintf(file_convergenza_edw, "%d;%.15E;%.15E\n", j, pot, rel_error);

            printf("{ED, Range:%d, RelErr:%.3E}: %.10E\n", j, rel_error, pot);
            if (rel_error == 0) break;
            last_pot = pot;
        }
    }

    //===== COULOMB =====//
    FILE *file_convergenza_coulomb = fopen("./data/convergenza_range/range_variabile_coulomb.csv", "w");
    last_pot = 0;
    printf("---------------------\n");
    for (size_t i = 0; i < _MAX_RANGE; i++)
    {
        _R_RANGE = i;
        double pot = getCoulombPotential(&system);
        double rel_error = fabs((last_pot - pot) / pot);
        fprintf(file_convergenza_coulomb, "%d;%.15E;%.15E\n", i, pot, rel_error);

        printf("{C, Range:%d, RelErr:%.3E}: %.10E\n", i, rel_error, pot);
        if (rel_error == 0) break;
        last_pot = pot;
    }
    free(system.particles);
}
