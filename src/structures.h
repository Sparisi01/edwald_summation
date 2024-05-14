#ifndef STRUCTURES_H
#define STRUCTURES_H

typedef struct Vec2
{
    double x, y;
} Vec2;

typedef struct Vec3
{
    double x, y, z;
} Vec3;

typedef struct Vec4
{
    double x, y, z, w;
} Vec4;

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

typedef struct PartialForStruct
{
    int start_index;
    int end_index;
    int n_particles;
    Particle *particles;
    Vec3 *tmp_forces;
} PartialForStruct;

#endif