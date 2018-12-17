#include <assert.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#include "data.h"
#include "io.h"
#include "traj.h"
#include "util.h"

static char buf[BUFSIZ];
static double *bufd;
static int bufdlen;

void init_buf(struct data *data)
{
    //id type mol mass x y z
    bufdlen = 1 + (data->nmols * ((6 * data->molsize) + 1));
    bufd = malloc(bufdlen * sizeof(double));
    assert(bufd != NULL);
}

void free_buf(void)
{
    free(bufd);
}

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
            data->xhi = hi;
            data->xlo = lo;
            data->xlen = fabs(hi - lo);
            continue;
        }
        if (strstr(buf, "ylo")) {
            double lo, hi;
            sscanf(buf, "%lf %lf", &lo, &hi);
            data->yhi = hi;
            data->ylo = lo;
            data->ylen = fabs(hi - lo);
            continue;
        }
        if (strstr(buf, "zlo")) {
            double lo, hi;
            sscanf(buf, "%lf %lf", &lo, &hi);
            data->zhi = hi;
            data->zlo = lo;
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
        frame->mol[mol].id = mol;

        R[id][0] = x;
        R[id][1] = y;
        R[id][2] = z;

        types[id] = type;
    }

    for (int i = 0; i < data->nmols; ++i) {
        assert(frame->mol[i].m == frame->mol[i].n);
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

long int read_frame_bin(gzFile fp, struct frame *frame, struct data *data)
{
    /* Read from file */
    int nread = gzread(fp, bufd, bufdlen * sizeof(double));
    if (nread < bufdlen) {
        if (gzeof(fp))
            return -1;
        else
            exit(EXIT_FAILURE);
    }

    /* Unpack buffer */
    double *b = bufd;
    long int step = (long int)(*b++);
    frame->step = step;
    for (int i = 0; i < data->nmols; ++i) {
        frame->mol[i].id = (int)(*b++);
        frame->mol[i].m = data->molsize;
        frame->mol[i].n = data->molsize;
        for (int j = 0; j < frame->mol[i].m; ++j) {
            frame->mol[i].atoms[j] = (int)(*b++);
            frame->mol[i].types[j] = (int)(*b++);
            frame->mol[i].mass[j] = *b++;
            frame->mol[i].R[j][0] = *b++;
            frame->mol[i].R[j][1] = *b++;
            frame->mol[i].R[j][2] = *b++;
        }
    }

    return step;
}

void write_frame_bin(gzFile fp, struct frame *frame, struct data *data)
{
    double *b = bufd;

    /* Build buffer */
    *b++ = (double)frame->step;
    for (int i = 0; i < data->nmols; ++i) {
        //mol id type mass x y z
        *b++ = (double)frame->mol[i].id;
        for (int j = 0; j < frame->mol[i].m; ++j) {
            *b++ = (double)frame->mol[i].atoms[j];
            *b++ = (double)frame->mol[i].types[j];
            *b++ = frame->mol[i].mass[j];
            *b++ = frame->mol[i].R[j][0];
            *b++ = frame->mol[i].R[j][1];
            *b++ = frame->mol[i].R[j][2];
        }
    }
    int len = b - bufd;
    assert(len == bufdlen);

    /* Write it */
    gzwrite(fp, bufd, len * sizeof(double));
}

void write_frame(gzFile fp, struct frame *frame, struct data *data)
{
    sprintf(buf, "ITEM: TIMESTEP\n");
    gzputs(fp, buf);

    sprintf(buf, "%ld\n", frame->step);
    gzputs(fp, buf);

    sprintf(buf, "ITEM: NUMBER OF ATOMS\n");
    gzputs(fp, buf);

    sprintf(buf, "%d\n", data->natoms);
    gzputs(fp, buf);

    sprintf(buf, "ITEM: BOX BOUNDS pp pp pp\n");
    gzputs(fp, buf);

    sprintf(buf, "%-16.10lf %-16.10lf\n", data->xlo, data->xhi);
    gzputs(fp, buf);

    sprintf(buf, "%-16.10lf %-16.10lf\n", data->ylo, data->yhi);
    gzputs(fp, buf);

    sprintf(buf, "%-16.10lf %-16.10lf\n", data->zlo, data->zhi);
    gzputs(fp, buf);

    sprintf(buf, "ITEM: ATOMS id type mol xu yu zu\n");
    gzputs(fp, buf);

    int id, type, mol;
    double x, y, z;
    for (int i = 0; i < data->nmols; ++i) {
        for (int j = 0; j < frame->mol[i].m; ++j) {
            id   = frame->mol[i].atoms[j] + 1;
            type = frame->mol[i].types[j];
            mol  = frame->mol[i].id + 1;
            x    = frame->mol[i].R[j][0];
            y    = frame->mol[i].R[j][1];
            z    = frame->mol[i].R[j][2];
            sprintf(buf, "%d %d %d %.6f %.6f %.6f\n", id, type, mol, x, y, z);
            gzputs(fp, buf);
        }
    }

}

void write_entropy(char *fn, struct data *data, double *S)
{
    FILE *fp = fopen(fn, "w");
    if (fp == NULL)
        fprintf(stderr, "error: %s\n", strerror(errno));

    printf("writing quasi-harmonic entropy\n");

    for (int i = 0; i < data->nmols; ++i) {
        if (fp)
            fprintf(fp, "%20.10f\n", S[i]);
        else
            printf("%20.10f\n", S[i]);
    }

    if (fp)
        fclose(fp);
}
