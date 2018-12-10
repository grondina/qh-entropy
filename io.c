#include <assert.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "io.h"
#include "data.h"

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

int read_frame(void)
{
    return 0;
}

void write_frame(void)
{
}
