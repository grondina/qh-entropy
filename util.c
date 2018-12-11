#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "data.h"
#include "util.h"

int cmpint(const void *p1, const void *p2)
{
    int a = *((int *)p1);
    int b = *((int *)p2);
    return a - b;
}

static double distance(double x1, double y1, double z1,
                       double x2, double y2, double z2)
{
    return sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2));
}

double gyration(struct molecule *mol)
{
    double sum = 0;
    double xi, yi, zi;
    double xj, yj, zj;

    for (int i = 0; i < mol->m; ++i) {

        xi = mol->R[i][0];
        yi = mol->R[i][1];
        zi = mol->R[i][2];

        for (int j = 0; j < mol->m; ++j) {

            xj = mol->R[j][0];
            yj = mol->R[j][1];
            zj = mol->R[j][2];

            double dist = distance(xi, yi, zi, xj, yj, zj);
            sum += dist*dist;
        }
    }

    return sum/(2.0 * mol->m * mol->m);
}
