#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define PI 3.1415926535897932384626433 // An over-accurate PI
#define SQR_2 1.4142135623730950488    // Square root of 2
#define SQR_3 1.73205080757            // Square root of 3
#define SQR_PI 1.77245385091           // Square root of PI
#define CELL_LENGHT 2.
#define CELL_VOLUME (CELL_LENGHT * CELL_LENGHT * CELL_LENGHT)
#define RAND_SEED 5
/* #define ALPHA (5.6 / CELL_LENGHT)
#define SIGMA (1 / (SQR_2 * ALPHA)) */
#define ALPHA 1e-6
#define N_K_RANGE 0
#define FORCE_TYPE_CONSTANT 1 // F = FORCE_TYPE_CONSTANT * (Q1*Q2)/r^2

#include "statistic.h"

struct Vec3
{
    double x, y, z;
};

struct Particle
{
    double x, y, z;
    double vx, vy, vz;
    double mass;
    double charge;
};

struct System
{
    double time;
    int n_particles;
    struct Particle *particles;
    struct Vec3 *forces;
};

struct Observables
{
    double *time;
    double *temperature;
    double *pressure;
    double *kinetic_energy;
    double *potential_energy;
};

int writeParticlesPositions(struct Particle *particles, int n_particles, FILE *file)
{
    if (!file) return 1;
    if (!particles) return 1;

    for (size_t i = 0; i < n_particles; i++)
    {
        fprintf(file, "%.5E %.5E %.5E\n", particles[i].x, particles[i].y, particles[i].z);
    }
    return 0;
}

double temperature(double *vel)
{
    return 1;
}

double kinetic_energy(struct Particle *particles, int n_particles)
{
    double kinetic_energy = 0;

    for (size_t i = 0; i < n_particles; i++)
    {
        kinetic_energy += 0.5 * particles[i].mass * (particles[i].vx * particles[i].vx);
        kinetic_energy += 0.5 * particles[i].mass * (particles[i].vy * particles[i].vy);
        kinetic_energy += 0.5 * particles[i].mass * (particles[i].vz * particles[i].vz);
    }

    return kinetic_energy;
}

double potential_energy(struct Particle *particles, int n_particles)
{
    double pot_energy = 0;
    for (size_t i = 0; i < n_particles - 1; i++)
    {
        for (size_t j = i + 1; j < n_particles; j++)
        {
            double r_ij_x = particles[i].x - particles[j].x - rint((particles[i].x - particles[j].x) / CELL_LENGHT) * CELL_LENGHT; //
            double r_ij_y = particles[i].y - particles[j].y - rint((particles[i].y - particles[j].y) / CELL_LENGHT) * CELL_LENGHT; //
            double r_ij_z = particles[i].z - particles[j].z - rint((particles[i].z - particles[j].z) / CELL_LENGHT) * CELL_LENGHT; //
            double r_ij = sqrt(r_ij_x * r_ij_x + r_ij_y * r_ij_y + r_ij_z * r_ij_z);

            pot_energy += FORCE_TYPE_CONSTANT * particles[i].charge * particles[j].charge / r_ij;
        }
    }
    return pot_energy;
}

int restore_positions_in_lattice_first_cell(struct Particle *particles, int n_particles)
{
    for (size_t i = 0; i < n_particles; i++)
    {
        particles[i].x = particles[i].x - rint(particles[i].x / CELL_LENGHT) * CELL_LENGHT;
        particles[i].y = particles[i].y - rint(particles[i].y / CELL_LENGHT) * CELL_LENGHT;
        particles[i].z = particles[i].z - rint(particles[i].z / CELL_LENGHT) * CELL_LENGHT;
    }
    return 0;
}

struct Vec3 *edwald_summation(struct System *system, double *args)
{
    struct Vec3 *tmp_forces = (struct Vec3 *)malloc(system->n_particles * sizeof(struct Vec3));
    if (!tmp_forces)
    {
        return NULL;
    }

    for (size_t i = 0; i < system->n_particles; i++)
    {
        tmp_forces[i].x = 0;
        tmp_forces[i].y = 0;
        tmp_forces[i].z = 0;
    }

    // Bring back all the particles in the (0,0,0) cell
    restore_positions_in_lattice_first_cell(system->particles, system->n_particles);

    const int SELECTED_PARTICLE = 0;
    static int r_k_print_countdown = 5;
    double real_sum = 0;
    double k_sum = 0;

    for (size_t i = 0; i < system->n_particles - 1; i++)
    {
        for (size_t j = i + 1; j < system->n_particles; j++)
        {
            double r_ij_x = system->particles[i].x - system->particles[j].x - rint((system->particles[i].x - system->particles[j].x) / CELL_LENGHT) * CELL_LENGHT; //
            double r_ij_y = system->particles[i].y - system->particles[j].y - rint((system->particles[i].y - system->particles[j].y) / CELL_LENGHT) * CELL_LENGHT; //
            double r_ij_z = system->particles[i].z - system->particles[j].z - rint((system->particles[i].z - system->particles[j].z) / CELL_LENGHT) * CELL_LENGHT; //
            double r_ij = sqrt(r_ij_x * r_ij_x + r_ij_y * r_ij_y + r_ij_z * r_ij_z);                                                                               // |r_i - r_j|

            /* if (r_ij >= CELL_L / 2) {
                continue;  // CUTOF potential
            } */

            double classical_force = system->particles[i].charge * system->particles[j].charge / (r_ij * r_ij * r_ij);         // Classic Coulomb Like force
            double edwald_correction = 2 * ALPHA / SQR_PI * exp(-(ALPHA * r_ij) * (ALPHA * r_ij)) * r_ij + erfc(ALPHA * r_ij); // Edwald correction in real space

            // Force on particle i due to j
            tmp_forces[i].x += FORCE_TYPE_CONSTANT * r_ij_x * classical_force * edwald_correction;
            tmp_forces[i].y += FORCE_TYPE_CONSTANT * r_ij_y * classical_force * edwald_correction;
            tmp_forces[i].z += FORCE_TYPE_CONSTANT * r_ij_z * classical_force * edwald_correction;

            // tmp_forces on particle j due to i using third principle
            tmp_forces[j].x -= FORCE_TYPE_CONSTANT * r_ij_x * classical_force * edwald_correction;
            tmp_forces[j].y -= FORCE_TYPE_CONSTANT * r_ij_y * classical_force * edwald_correction;
            tmp_forces[j].z -= FORCE_TYPE_CONSTANT * r_ij_z * classical_force * edwald_correction;

            if (i == SELECTED_PARTICLE)
            {
                real_sum += FORCE_TYPE_CONSTANT * r_ij_x * classical_force * edwald_correction;
                real_sum += FORCE_TYPE_CONSTANT * r_ij_y * classical_force * edwald_correction;
                real_sum += FORCE_TYPE_CONSTANT * r_ij_z * classical_force * edwald_correction;
            }
        }

        // SECTION - Reciprocal space force

        /**NOTE - La formula non prevede l'esclusione dell'interazione tra particella i e j
         * se esse appartengono a celle diverse. Noto però che in questo caso il sin da valore
         * nullo. Posso escludere l'interazione.
         **/
        for (size_t j = i + 1; j < system->n_particles; j++)
        {
            for (int n_k_x = -N_K_RANGE; n_k_x <= N_K_RANGE; n_k_x++)
            {
                for (int n_k_y = -N_K_RANGE; n_k_y <= N_K_RANGE; n_k_y++)
                {
                    for (int n_k_z = -N_K_RANGE; n_k_z <= N_K_RANGE; n_k_z++)
                    {
                        if (n_k_x == 0 && n_k_y == 0 && n_k_z == 0)
                        {
                            continue;
                        }
                        double k_x = n_k_x * 2 * PI / CELL_LENGHT;               // Componente x vettore spazione reciproco
                        double k_y = n_k_y * 2 * PI / CELL_LENGHT;               // Componente y vettore spazione reciproco
                        double k_z = n_k_z * 2 * PI / CELL_LENGHT;               // Componente z vettore spazione reciproco
                        double k_vec_mag2 = (k_x * k_x + k_y * k_y + k_z * k_z); // Magnitude square of reciprocal-lattice vector

                        double dot_prod_pos_k = (k_x * (system->particles[i].x - system->particles[j].x) + k_y * (system->particles[i].y - system->particles[j].y) + k_z * (system->particles[i].z - system->particles[j].z));
                        double long_distance_force = 4 * PI / CELL_VOLUME * system->particles[i].charge * system->particles[j].charge;
                        long_distance_force *= exp(-k_vec_mag2 / (2 * ALPHA * ALPHA)) / k_vec_mag2 * sin(dot_prod_pos_k);

                        // tmp_forces on particle i due to j copies
                        tmp_forces[i].x += FORCE_TYPE_CONSTANT * k_x * long_distance_force;
                        tmp_forces[i].z += FORCE_TYPE_CONSTANT * k_y * long_distance_force;
                        tmp_forces[i].y += FORCE_TYPE_CONSTANT * k_z * long_distance_force;
                        // tmp_forces on particle j due to i copies using third principle
                        tmp_forces[j].x -= FORCE_TYPE_CONSTANT * k_x * long_distance_force;
                        tmp_forces[j].y -= FORCE_TYPE_CONSTANT * k_y * long_distance_force;
                        tmp_forces[j].z -= FORCE_TYPE_CONSTANT * k_z * long_distance_force;

                        if (i == SELECTED_PARTICLE)
                        {
                            k_sum += FORCE_TYPE_CONSTANT * k_x * long_distance_force;
                            k_sum += FORCE_TYPE_CONSTANT * k_y * long_distance_force;
                            k_sum += FORCE_TYPE_CONSTANT * k_z * long_distance_force;
                        }
                    }
                }
            }
        }
        /*!SECTION */
    }
    if (r_k_print_countdown > 0)
    {
        printf("REAL SUM: %.10E\nRECIPROCAL SUM: %.10E\n", real_sum, k_sum);
        r_k_print_countdown--;
    }

    return tmp_forces;
}

int verletPropagationStep(struct System *system, double time_step, struct Vec3 *(*forceFunction)(struct System *system, double *args), double *args)
{

    for (size_t i = 0; i < system->n_particles; i++)
    {
        system->particles[i].x += system->particles[i].vx * time_step + 0.5 / system->particles[i].mass * system->forces[i].x * time_step * time_step;
        system->particles[i].y += system->particles[i].vy * time_step + 0.5 / system->particles[i].mass * system->forces[i].y * time_step * time_step;
        system->particles[i].z += system->particles[i].vz * time_step + 0.5 / system->particles[i].mass * system->forces[i].z * time_step * time_step;
    }

    system->time += time_step;

    struct Vec3 *new_forces = forceFunction(system, args);
    if (!new_forces)
    {
        return 1;
    }

    for (size_t i = 0; i < system->n_particles; i++)
    {
        system->particles[i].vx += 0.5 / system->particles[i].mass * (system->forces[i].x + new_forces[i].x) * time_step;
        system->particles[i].vy += 0.5 / system->particles[i].mass * (system->forces[i].y + new_forces[i].y) * time_step;
        system->particles[i].vz += 0.5 / system->particles[i].mass * (system->forces[i].z + new_forces[i].z) * time_step;
    }

    free(system->forces);
    system->forces = new_forces;

    return 0;
}

int main(int argc, char const *argv[])
{

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
    system.n_particles = 4;
    system.particles = (struct Particle *)malloc(sizeof(struct Particle) * system.n_particles);
    if (!system.particles)
    {
        perror("Error, malloc 'system.particles' returned NULL:");
        exit(EXIT_FAILURE);
    }

    //----------------------------------------------
    double start_time = 0;
    double end_time = 30;
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
        // system.particles[i].x = rand() / (RAND_MAX + 1.0) * CELL_LENGHT - CELL_LENGHT / 2;
        // system.particles[i].y = rand() / (RAND_MAX + 1.0) * CELL_LENGHT - CELL_LENGHT / 2;
        // system.particles[i].z = rand() / (RAND_MAX + 1.0) * CELL_LENGHT - CELL_LENGHT / 2;

        system.particles[i].x = 0;
        system.particles[i].y = 0;
        system.particles[i].z = 0;

        system.particles[i].vx = 0;
        system.particles[i].vy = 0;
        system.particles[i].vz = 0;

        system.particles[i].mass = 1;
        system.particles[i].charge = 0.01;
    }

    system.particles[0].x = 0.3;
    system.particles[1].x = -0.3;
    system.particles[2].y = 0.3;
    system.particles[3].y = -0.3;

    system.particles[0].charge = 0.01;
    system.particles[1].charge = 0.01;
    system.particles[2].charge = 0.01;
    system.particles[3].charge = 0.01;

    system.particles[0].mass = 1;
    system.particles[1].mass = 1;
    system.particles[2].mass = 1;
    system.particles[3].mass = 1;

    system.forces = edwald_summation(&system, NULL);

    // -----------------------------

    printf("--------------------\n");
    printf("AVVIO SIMULAZIONE\n");
    printf("--------------------\n");
    printf("N Particelle: %d\n", system.n_particles);
    printf("Volume: %lf\n", CELL_VOLUME);
    printf("Δt, dt: %lf, %lf\n", end_time - start_time, time_step);
    printf("--------------------\n");

    writeParticlesPositions(system.particles, system.n_particles, file_start_particles_pos);

    int n_time_measure = 50;
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
        observables.pressure = 0;
        observables.temperature = 0;

        if (i < n_time_measure)
        {
            sum_comput_times += clock() - comput_time_start;
        }
        if (i == n_time_measure - 1)
        {
            double time_spent = (double)sum_comput_times / n_time_measure;
            int tot_time_sec = time_spent * n_time_step / CLOCKS_PER_SEC;
            int sec = tot_time_sec % 60;        // Estimated time seconds
            int min = (tot_time_sec / 60) % 60; // Estimated time minutes
            int h = tot_time_sec / 3600;        // Estimated time hours
            printf("--------------------\nVerlet time: %d ms\nTempo totale stimato: %d:%d:%d [h:m:s]\n--------------------\n", (int)time_spent, h, min, sec);
        }
    }

    writeParticlesPositions(system.particles, system.n_particles, file_end_particles_pos);

    // Save energies in file and print statistics
    for (size_t i = 0; i < n_time_step; i++)
    {
        fprintf(file_thermodinamic, "%.10E %.10E %.10E %.10E\n", start_time + i * time_step, observables.kinetic_energy[i] + observables.potential_energy[i], observables.kinetic_energy[i], observables.potential_energy[i]);
    }

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
