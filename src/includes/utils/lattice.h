#if !defined(LATTICE_H)
#define LATTICE_H

#include "../constants.h"
#include "../structures.h"
#include <math.h>

typedef struct LatticeCell
{
    int n_particle;
    Particle *particles;
    Vec3 *lattice_vectors;
} LatticeCell;

int initialize_in_lattice_by_density(System *s, LatticeCell cell, int n_cells, double density)
{
    s->n_particles = cell.n_particle * pow(n_cells, 3);
    s->cell_lenght = cbrt(s->n_particles / density);
    double lattice_cell_lenght = s->cell_lenght / n_cells;

    if (s->particles)
    {
        free(s->particles);
    }

    s->particles = (Particle *)malloc(sizeof(Particle) * s->n_particles);
    if (!s->particles)
    {
        return 1;
    }

    return 0;
}

int restore_positions_in_lattice_first_cell(Particle *particles, int n_particles, double cell_lenght)
{
    for (size_t i = 0; i < n_particles; i++)
    {
        particles[i].x = particles[i].x - rint(particles[i].x / cell_lenght) * cell_lenght;
        particles[i].y = particles[i].y - rint(particles[i].y / cell_lenght) * cell_lenght;
        particles[i].z = particles[i].z - rint(particles[i].z / cell_lenght) * cell_lenght;
    }
    return 0;
}

/// @brief Apply the minimum image convention.
/// Link: https://en.wikipedia.org/wiki/Periodic_boundary_conditions
/// @param pos
/// @param cell_length
double minimum_image(double pos, double cell_length)
{
    return pos - rint(pos / cell_length) * cell_length;
}

#endif // LATTICE_UTILS_H
