#if !defined(LATTICE_CELLS_H)
#define LATTICE_CELLS_H

#include "utils/lattice.h"

LatticeCell NaOh = {
    .n_particle = {},
    .particles = {},
    .lattice_vectors = {},
};

LatticeCell FCC = {
    .n_particle = 4,
    .particles = {},
    .lattice_vectors = {},
};

LatticeCell BCC = {
    .n_particle = 2,
    .particles = {},
    .lattice_vectors = {},
};

#endif // LATTICE_CELLS_H
