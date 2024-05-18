#include "structures.h"

#ifndef INTERPOLATION_H
#define INTERPOLATION_H

double lerp1D(
    double x,
    double x1,
    double y1,
    double x2,
    double y2)
{
    return y1 + (x - x1) * (y2 - y1) * (x2 - x1) / (x2 * x2 - x1 * x1);
}

double lerp2D(
    double x,
    double y,
    Vec3 vert1, Vec3 vert2,
    Vec3 vert3, Vec3 vert4)
{
    double l1 = lerp1D(x, vert1.x, vert1.z, vert2.x, vert2.z);
    double l3 = lerp1D(x, vert3.x, vert3.z, vert4.x, vert4.z);
    return lerp1D(y, vert1.y, l1, vert3.y, l3);
}

#define TOVEC3(vec4)           \
    (Vec3)                     \
    {                          \
        vec4.x, vec4.y, vec4.w \
    }

// Formula for trilinear interpolation.
// Source https://en.wikipedia.org/wiki/Trilinear_interpolation#/media/File:3D_interpolation2.svg
// Vertex must be in the current position (x,y,z):
// 1 - (0,0,0) // 4 - (1,1,0) // 7 - (1,0,1)
// 2 - (1,0,0) // 5 - (0,0,1) // 8 - (1,1,1)
// 3 - (0,1,0) // 6 - (0,1,1)
double lerp3D(
    double x,
    double y,
    double z,
    Vec4 vert1, Vec4 vert2, Vec4 vert3, Vec4 vert4,
    Vec4 vert5, Vec4 vert6, Vec4 vert7, Vec4 vert8)
{
    // Vertex Plane 1 (1,2,3,4)
    double l1 = lerp2D(x, y, TOVEC3(vert1), TOVEC3(vert2), TOVEC3(vert3), TOVEC3(vert4));

    // Vertex Plane 2 (5,6,7,8)
    double l2 = lerp2D(x, y, TOVEC3(vert5), TOVEC3(vert6), TOVEC3(vert7), TOVEC3(vert8));

    // Linear interpolation between the last two
    return lerp1D(z, vert1.z, l1, vert5.z, l2);
}

#endif