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
#include "util.h"

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

    /* Reset read counter for all molecules */
    for (int i = 0; i < frame->nmols; ++i)
        frame->mol[i].m = 0;

    /* TIMESTEP */
    if (gzgets(fp, buf, BUFSIZ) == NULL) {
        if (gzeof(fp)) {
            return -1;
        } else {
            int err;
            fprintf(stderr, "error: %s\n", gzerror(fp, &err));
            exit(EXIT_FAILURE);
        }
    }

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
    assert(gzgets(fp, buf, BUFSIZ) != NULL);

    /* Read atoms */

    double (*R)[3] = malloc(sizeof (double[data->natoms][3]));
    assert(R != NULL);

    int *types = malloc(data->natoms * sizeof(int));
    assert(types != NULL);

    for (int i = 0; i < data->natoms; ++i) {

        int id, type, mol;
        double x, y, z;

        assert(gzgets(fp, buf, BUFSIZ) != NULL);
        sscanf(buf, "%d %d %d %lf %lf %lf", &id, &type, &mol, &x, &y, &z);

        assert(id <= data->natoms);
        assert(mol <= data->nmols);

        mol--;
        id--;

        int j = frame->mol[mol].m;
        assert(j < frame->mol[mol].n);

        frame->mol[mol].m++;
        frame->mol[mol].atoms[j] = id;

        R[id][0] = x;
        R[id][1] = y;
        R[id][2] = z;

        types[id] = type;
    }

    for (int i = 0; i < data->nmols; ++i) {
        qsort(frame->mol[i].atoms, frame->mol[i].m, sizeof(int), cmpint);
        for (int j = 0; j < frame->mol[i].m; ++j) {
            int id = frame->mol[i].atoms[j];
            frame->mol[i].R[j][0] = R[id][0];
            frame->mol[i].R[j][1] = R[id][1];
            frame->mol[i].R[j][2] = R[id][2];
            frame->mol[i].types[j] = types[id];
            frame->mol[i].mass[j] = get_mass(data, types[id]);
        }
    }

    free(R);
    free(types);

    return step;
}

void write_frame(void)
{
}
