#ifndef STRUCTURES_H
#define STRUCTURES_H

typedef struct Vec3
{
    double x, y, z;
} Vec3;

typedef struct Particle
{
    double x, y, z;
    double vx, vy, vz;
    double mass;
    double charge;
} Particle;

typedef struct System
{
    double time;
    int n_particles;
    struct Particle *particles;
    struct Vec3 *forces;
} System;

typedef struct Observables
{
    double *time;
    double *temperature;
    double *pressure;
    double *kinetic_energy;
    double *potential_energy;
} Observables;

#endif