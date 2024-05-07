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

/* double lerp3D(
    Vec3 x, Vec3 x1, Vec3 x2, Vec3 x3, Vec3 x4, Vec3 x5, Vec3 x6, Vec3 x7, Vec3 x8,
    double v1, double v2, double v3, double v4, double v5, double v6, double v7, double v8)
{

    double w1 = (x.x - x1.x) * (x.y - x1.y) * (x.z - x1.z) * v8;
    double w2 = (x2.x - x.x) * (x.y - x2.y) * (x.z - x2.z) * v7;
    double w3 = (x.x - x3.x) * (x3.y - x.y) * (x.z - x3.z) * v6;
    double w4 = (x4.x - x.x) * (x4.y - x.y) * (x4.z - x.z) * v5;
    double w5 = (x.x - x5.x) * (x.y - x5.y) * (x5.z - x.z) * v4;
    double w6 = (x6.x - x.x) * (x.y - x6.y) * (x6.z - x.z) * v3;
    double w7 = (x.x - x7.x) * (x7.y - x.y) * (x7.z - x.z) * v2;
    double w8 = (x8.x - x.x) * (x8.y - x.y) * (x8.z - x.z) * v1;

    return (w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8) / ((x8.x - x1.x) * (x8.y - x1.y) * (x8.z - x1.z));
} */

#endif