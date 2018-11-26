#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <string.h>
#include <errno.h>

/* 2D array indexing for row major order */
#define ARRAY2D(A, n1, n2, N2) ((A)[(n2) + ((N2)*(n1))])

#define NCHAINS  500
#define CHAINLEN 8

/* chain */
typedef struct {
    int n;       /* total number of atoms in the chain */
    int m;       /* number read so far */
    double *r;   /* coordinate matrix (n x 3), row major order  */
    double *x;
    double *y;
    double *z;
} chain_t;

char buf[BUFSIZ];

char *read_line(gzFile fp)
{
    char *s;

    s = gzgets(fp, buf, BUFSIZ);
    if (s == NULL) {
        fprintf(stderr, "Error: could not read from file: %s\n", strerror(errno));
        exit(1);
    }

    return s;
}

int read_frame(gzFile fp, chain_t *chains, int nchains)
{
    int i, j;
    int step, natoms;
    int id, mol;
    double x, y, z;

    /* TIMESTEP */
    if (gzgets(fp, buf, BUFSIZ) == NULL) {
        return -1;
    }

    /* TIMESTEP value */
    read_line(fp);
    sscanf(buf, "%d", &step);

    /* NUMBER OF ATOMS */
    read_line(fp);

    /* NUMBER OF ATOMS value */
    read_line(fp);
    sscanf(buf, "%d", &natoms);

    /* BOX BOUNDS pp pp pp */
    read_line(fp);
    read_line(fp);
    read_line(fp);
    read_line(fp);

    /* ATOMS id mol xu yu zu */
    read_line(fp);

    /* Read atoms */
    for (i = 0; i < natoms; i++) {
        read_line(fp);
        sscanf(buf, "%d %d %lf %lf %lf", &id, &mol, &x, &y, &z);
        assert(id <= natoms);
        assert(mol <= nchains);
        mol--;
        j = chains[mol].m;
        chains[mol].x[j] = x;
        chains[mol].y[j] = y;
        chains[mol].z[j] = z;
        chains[mol].m++;
        assert(chains[mol].m <= chains[mol].n);
    }

    for (i = 0; i < nchains; i++) {
        chains[i].m = 0;
    }

    return step;
}

double read_mass(char *fn)
{
    int type;
    double mass;
    FILE *fp;

    if (fn == NULL) {
        return -1;
    }

    fp = fopen(fn, "r");
    if (fp == NULL) {
        fprintf(stderr, "Error: could not open file %s\n", fn);
        exit(1);
    }

    while (fgets(buf, BUFSIZ, fp)) {
        if (strstr(buf, "Masses")) {
            assert(fgets(buf, BUFSIZ, fp));
            assert(fgets(buf, BUFSIZ, fp));
            sscanf(buf, "%d %lf", &type, &mass);
            fclose(fp);
            return mass;
        }
    }

    fclose(fp);
    return -1;

}

int main(int argc, char *argv[])
{
    int i, j;
    int step;
    chain_t chains[NCHAINS];
    double clen = CHAINLEN;
    double mass;
    gzFile fp;
    double xcom, ycom, zcom;

    for (i = 0; i < NCHAINS; i++) {
        chains[i].n = clen;
        chains[i].m = 0;
        chains[i].r = NULL;
        chains[i].x = malloc(clen * sizeof(double));
        chains[i].y = malloc(clen * sizeof(double));
        chains[i].z = malloc(clen * sizeof(double));
        assert(chains[i].x != NULL);
        assert(chains[i].y != NULL);
        assert(chains[i].z != NULL);
    }

    fp = gzopen("traj.dump.gz", "r");
    assert(fp != NULL);

    if (gzbuffer(fp, 0xF0000) == -1) {
        fprintf(stderr, "error: gzbuffer\n");
        return EXIT_FAILURE;
    }

    mass = read_mass("data.lmp");
    assert(mass > 0);

    while ((step = read_frame(fp, chains, NCHAINS)) >= 0) {

        printf("%d\n", step);

        /* Remove COM from each molecule */
        for (i = 0; i < NCHAINS; i++) {
            xcom = 0;
            ycom = 0;
            zcom = 0;
            for (j = 0; j < chains[i].n; j++) {
                xcom += (chains[i].x[j] * mass);
                ycom += (chains[i].y[j] * mass);
                zcom += (chains[i].z[j] * mass);
            }
            xcom /= (chains[i].n * mass);
            ycom /= (chains[i].n * mass);
            zcom /= (chains[i].n * mass);
            for (j = 0; j < chains[i].n; j++) {
                chains[i].x[j] -= xcom;
                chains[i].y[j] -= ycom;
                chains[i].z[j] -= zcom;
            }
        }
        
        break;
    }

    gzclose(fp);

    return 0;
}
