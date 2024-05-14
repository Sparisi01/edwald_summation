#include "interpolation.h"
#include "structures.h"

#include <stdio.h>

int main(int argc, char const *argv[])
{
    // Linear interpolation test
    double x1 = 0.;
    double y1 = 0.;
    double x2 = 1.;
    double y2 = 1.;

    printf("TEST lerp1D %lf\n", lerp1D(1, x1, y1, x2, y2));
    printf("TEST lerp1D %lf\n", lerp1D(0.5, x1, y1, x2, y2));
    printf("TEST lerp1D %lf\n", lerp1D(0, x1, y1, x2, y2));

    // Bilinear interpolation
    Vec3 vertex1 = {0, 0, 0};
    Vec3 vertex2 = {1, 0, 0.5};
    Vec3 vertex3 = {0, 1, 0.5};
    Vec3 vertex4 = {1, 1, 1};

    printf("TEST lerp2D %lf\n", lerp2D(0, 0, vertex1, vertex2, vertex3, vertex4));
    printf("TEST lerp2D %lf\n", lerp2D(0.5, 0.5, vertex1, vertex2, vertex3, vertex4));
    printf("TEST lerp2D %lf\n", lerp2D(1, 1, vertex1, vertex2, vertex3, vertex4));

    // Trilinear interpolation

    Vec4 v1 = {0, 0, 0, 0};
    Vec4 v2 = {1, 0, 0, 0};
    Vec4 v3 = {0, 1, 0, 0};
    Vec4 v4 = {1, 1, 0, 0};
    Vec4 v5 = {0, 0, 1, 6};
    Vec4 v6 = {1, 0, 1, 6};
    Vec4 v7 = {0, 1, 1, 6};
    Vec4 v8 = {1, 1, 1, 6};

    // Test Vertecies
    printf("TEST lerp3D %lf\n", lerp3D(0, 0, 0, v1, v2, v3, v4, v5, v6, v7, v8));
    printf("TEST lerp3D %lf\n", lerp3D(1, 1, 1, v1, v2, v3, v4, v5, v6, v7, v8));
    // In the current setup must be independent on x,y plane
    printf("TEST lerp3D %lf\n", lerp3D(0.5, 0.5, 0.5, v1, v2, v3, v4, v5, v6, v7, v8));
    printf("TEST lerp3D %lf\n", lerp3D(0.6, 0.6, 0.6, v1, v2, v3, v4, v5, v6, v7, v8));
    printf("TEST lerp3D %lf\n", lerp3D(0.6, 0.6, 0.7, v1, v2, v3, v4, v5, v6, v7, v8));

    return 0;
}