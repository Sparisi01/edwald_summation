
int verletPropagationStep(double *pos, double *vel, double **forces_ptr, double *masses, int n_particles, double t, double dt, double *(*F)(double t, double *pos, double *vel, int n_particles, void *args), void *args) {
    // Aggiorna posizioni
    double *forces = *forces_ptr;
    for (size_t i = 0; i < n_particles * 3; i += 3) {
        pos[i + 0] += vel[i + 0] * dt + 0.5 / masses[i / 3] * forces[i + 0] * dt * dt;
        pos[i + 1] += vel[i + 1] * dt + 0.5 / masses[i / 3] * forces[i + 1] * dt * dt;
        pos[i + 2] += vel[i + 2] * dt + 0.5 / masses[i / 3] * forces[i + 2] * dt * dt;
    }

    // Calcola forze agenti sulle particelle
    double *new_forces = F(t + dt, pos, vel, n_particles, args);
    if (!new_forces) {
        return 1;
    }

    // Aggiorna velocità
    for (size_t i = 0; i < n_particles * 3; i += 3) {
        vel[i + 0] += 0.5 / masses[i / 3] * (forces[i + 0] + new_forces[i + 0]) * dt;
        vel[i + 1] += 0.5 / masses[i / 3] * (forces[i + 1] + new_forces[i + 1]) * dt;
        vel[i + 2] += 0.5 / masses[i / 3] * (forces[i + 2] + new_forces[i + 2]) * dt;
    }

    // Sovrascrivi forze e libera memoria array vecchie forze
    free(*forces_ptr);
    *forces_ptr = new_forces;

    return 0;
}