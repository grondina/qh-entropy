#include <assert.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include "io.h"
#include "data.h"
#include "traj.h"

static char buf[BUFSIZ];

int read_data(char *fndata, struct data *data)
{
    assert(fndata != NULL);
    assert(data != NULL);

    FILE *fp;
    fp = fopen(fndata, "r");
    if (fp == NULL) {
        perror("read_data");
        return -1;
    }

    while (fgets(buf, BUFSIZ, fp)) {

        /* Number of atom types */
        if (strstr(buf, "atom types")) {
            int ntypes;
            sscanf(buf, "%d", &ntypes);
            data->ntypes = ntypes;
            continue;
        }

        /* Number of atoms */
        if (strstr(buf, "atoms")) {
            int natoms;
            sscanf(buf, "%d", &natoms);
            data->natoms = natoms;
            continue;
        }

        /* Atom masses */
        if (strstr(buf, "Masses")) {
            /* Make sure we know how many types are there */
            assert(data->ntypes > 0);
            /* Allocate vectors */
            data->mass = malloc(data->ntypes * sizeof(double));
            data->type = malloc(data->ntypes * sizeof(double));
            assert(data->mass != NULL);
            assert(data->type != NULL);
            /* Skip blank line */
            assert(fgets(buf, BUFSIZ, fp));
            /* Read ntypes masses */
            for (int i = 0; i < data->ntypes; ++i) {
                int type;
                double mass;
                assert(fgets(buf, BUFSIZ, fp));
                if (sscanf(buf, "%d %lf", &type, &mass) == 2) {
                    assert(type > 0);
                    data->type[i] = type;
                    data->mass[i] = mass;
                } else {
                    fprintf(stderr, "read_data: not enough masses\n");
                    return -1;
                }
            }
            continue;
        }

        /* Box dimensions */
        if (strstr(buf, "xlo")) {
            double lo, hi;
            sscanf(buf, "%lf %lf", &lo, &hi);
            data->xlen = fabs(hi - lo);
            continue;
        }
        if (strstr(buf, "ylo")) {
            double lo, hi;
            sscanf(buf, "%lf %lf", &lo, &hi);
            data->ylen = fabs(hi - lo);
            continue;
        }
        if (strstr(buf, "zlo")) {
            double lo, hi;
            sscanf(buf, "%lf %lf", &lo, &hi);
            data->zlen = fabs(hi - lo);
            continue;
        }
    }

    if (ferror(fp)) {
        perror("read_data");
        return -1;
    }

    fclose(fp);
    return 0;
}

long int read_frame(gzFile fp, struct frame *frame, struct data *data)
{
    assert(frame != NULL);
    assert(data != NULL);

    /* TIMESTEP */
    assert(gzgets(fp, buf, BUFSIZ) != NULL);

    /* TIMESTEP value */
    long int step;
    assert(gzgets(fp, buf, BUFSIZ) != NULL);
    sscanf(buf, "%ld", &step);
    frame->step = step;

    /* NUMBER OF ATOMS */
    assert(gzgets(fp, buf, BUFSIZ) != NULL);

    /* NUMBER OF ATOMS value */
    int natoms;
    assert(gzgets(fp, buf, BUFSIZ) != NULL);
    sscanf(buf, "%d", &natoms);
    assert(natoms == data->natoms);

    /* BOX BOUNDS pp pp pp */
    assert(gzgets(fp, buf, BUFSIZ) != NULL);

    /* x bounds */
    double lo, hi;
    assert(gzgets(fp, buf, BUFSIZ) != NULL);
    sscanf(buf, "%lf %lf", &lo, &hi);
    double xlen = fabs(hi - lo);

    /* y bounds */
    assert(gzgets(fp, buf, BUFSIZ) != NULL);
    sscanf(buf, "%lf %lf", &lo, &hi);
    double ylen = fabs(hi - lo);

    /* z bounds */
    assert(gzgets(fp, buf, BUFSIZ) != NULL);
    sscanf(buf, "%lf %lf", &lo, &hi);
    double zlen = fabs(hi - lo);

    /* Make sure box length isn't changing */
    assert(fabs(data->xlen - xlen) < 1.0e6);
    assert(fabs(data->ylen - ylen) < 1.0e6);
    assert(fabs(data->zlen - zlen) < 1.0e6);

    /* ATOMS id mol xu yu zu */
    assert(fabs(data->zlen - zlen) < 1.0e6);

    /* Read atoms */
    for (int i = 0; i < data->natoms; ++i) {
        int id, mol;
        double x, y, z;
        assert(gzgets(fp, buf, BUFSIZ) != NULL);
        sscanf(buf, "%d %d %lf %lf %lf", &id, &mol, &x, &y, &z);
        assert(id <= data->natoms);
        assert(mol <= data->nmols);
        mol--;
        int j = frame->mol[mol].m;
        frame->mol[mol].R[j][0] = x;
        frame->mol[mol].R[j][1] = y;
        frame->mol[mol].R[j][2] = z;
        frame->mol[mol].m++;
        assert(frame->mol[mol].m <= frame->mol[mol].m);
    }

    return step;
}

void write_frame(void)
{
}
