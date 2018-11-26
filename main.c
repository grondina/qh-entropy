#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <float.h>

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

inline double distance(double x1, double y1, double z1, double x2, double y2, double z2)
{
    return sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2));
}

inline double distance_sq(double x1, double y1, double z1, double x2, double y2, double z2)
{
    return ((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2));
}

void copy_chain(chain_t *c1, chain_t *c2)
{
    int i, j;

    assert(c1 != NULL);
    assert(c2 != NULL);

    c2->n = c1->n;
    c2->m = c1->m;

    if (c1->r != NULL) {
        if (c2->r == NULL) {
            c2->r = malloc(c2->n * 3 * sizeof(double));
            assert(c2->r != NULL);
        }
        for (i = 0; i < c2->n; i++) {
            for (j = 0; j < 3; j++) {
                ARRAY2D(c2->r, i, j, 3) = ARRAY2D(c1->r, i, j, 3);
            }
        }
    }

    if (c1->x != NULL) {
        if (c2->x == NULL) {
            c2->x = malloc(c2->n * sizeof(double));
            assert(c2->x != NULL);
        }
        for (i = 0; i < c2->n; i++) {
            c2->x[i] = c1->x[i];
        }
    }

    if (c1->y != NULL) {
        if (c2->y == NULL) {
            c2->y = malloc(c2->n * sizeof(double));
            assert(c2->y != NULL);
        }
        for (i = 0; i < c2->n; i++) {
            c2->y[i] = c1->y[i];
        }
    }

    if (c1->z != NULL) {
        if (c2->z == NULL) {
            c2->z = malloc(c2->n * sizeof(double));
            assert(c2->z != NULL);
        }
        for (i = 0; i < c2->n; i++) {
            c2->z[i] = c1->z[i];
        }
    }
}

void chain_convert(chain_t *c)
{
    int i;

    for (i = 0; i < c->n; i++) {
        ARRAY2D(c->r, i, 0, 3) = c->x[i];
        ARRAY2D(c->r, i, 1, 3) = c->y[i];
        ARRAY2D(c->r, i, 2, 3) = c->z[i];
    }
}

void chain_print_matrix(chain_t *c)
{
    int i, j;

    for (i = 0; i < c->n; i++) {
        for (j = 0; j < 3; j++) {
            printf("  % -6.4f%s", ARRAY2D(c->r, i, j, 3), (j == 2 ? "\n" : "  "));
        }
    }
}

int main(int argc, char *argv[])
{
    int i, j, k;
    int step;
    chain_t chains[NCHAINS];
    chain_t ref;
    double clen = CHAINLEN;
    double mass;
    gzFile fp;
    double xcom, ycom, zcom;
    double rg;
    double rg_min;
    double xj, yj, zj;
    double xk, yk, zk;

    for (i = 0; i < NCHAINS; i++) {
        chains[i].n = clen;
        chains[i].m = 0;
        chains[i].r = malloc(clen * 3 * sizeof(double));
        chains[i].x = malloc(clen * sizeof(double));
        chains[i].y = malloc(clen * sizeof(double));
        chains[i].z = malloc(clen * sizeof(double));
        assert(chains[i].r != NULL);
        assert(chains[i].x != NULL);
        assert(chains[i].y != NULL);
        assert(chains[i].z != NULL);
    }

    /* Reference */
    ref.n = clen;
    ref.m = 0;
    ref.r = malloc(clen * 3 * sizeof(double));
    ref.x = malloc(clen * sizeof(double));
    ref.y = malloc(clen * sizeof(double));
    ref.z = malloc(clen * sizeof(double));
    assert(ref.r != NULL);
    assert(ref.x != NULL);
    assert(ref.y != NULL);
    assert(ref.z != NULL);

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

        /* Save reference */
        rg_min = DBL_MAX;
        for (i = 0; i < NCHAINS; i++) {
            rg = 0;
            for (j = 0; j < chains[i].n; j++) {
                xj = chains[i].x[j];
                yj = chains[i].y[j];
                zj = chains[i].z[j];
                for (k = 0; k < chains[i].n; k++) {
                    xk = chains[i].x[k];
                    yk = chains[i].y[k];
                    zk = chains[i].z[k];
                    rg += distance_sq(xj, yj, zj, xk, yk, zk);
                }
                rg /= (2.0 * chains[i].n * chains[i].n);
                rg = sqrt(rg);
            }
            if (rg < rg_min) {
                rg_min = rg;
                copy_chain(&chains[i], &ref);
            }
        }

        /* Convert to matrix */
        for (i = 0; i < NCHAINS; i++) {
            chain_convert(&chains[i]);
            //printf("--- %3d ---\n", i);
            //chain_print_matrix(&chains[i]);
        }
        chain_convert(&ref);
        //printf("--- ref ---\n");
        //chain_print_matrix(&ref);
        
        break;
    }

    gzclose(fp);

    return 0;
}
