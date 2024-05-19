#if !defined(LATTICE_UTILS_H)
#define LATTICE_UTILS_H

#include "constants.h"
#include "structures.h"
#include <math.h>

typedef struct LatticeCell
{
    double length;
    int n_particle;
    Particle *particles;
    Vec3 *lattice_vectors;
} LatticeCell;

int initialize_in_lattice_by_cell(System *s, LatticeCell cell, int n_cells)
{
    s->n_particles = cell.n_particle * pow(n_cells, 3);
    s->cell_lenght = cell.length * n_cells;

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

#endif // LATTICE_UTILS_H
