#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <float.h>
#include <lapacke.h>
#include <lapacke_utils.h>
#include <cblas.h>

/* 2D array indexing for row major order */
#define ARRAY2D(A, n1, n2, N2) ((A)[(n2) + ((N2)*(n1))])

#define NCHAINS  500
#define CHAINLEN 8

/* chain */
typedef struct {
    int n;       /* total number of atoms in the chain */
    int m;       /* number read so far */
    double *r;   /* coordinate matrix (n x 3), row major order  */
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

void write_frame(gzFile fpo, chain_t *chains, int step, int nchains, int natoms, double *boxlo, double *boxhi)
{
    int i, j, id;
    double x, y, z;


    sprintf(buf, "ITEM: TIMESTEP\n");
    gzputs(fpo, buf);

    sprintf(buf, "%d\n", step);
    gzputs(fpo, buf);

    sprintf(buf, "ITEM: NUMBER OF ATOMS\n");
    gzputs(fpo, buf);

    sprintf(buf, "%d\n", natoms);
    gzputs(fpo, buf);

    sprintf(buf, "ITEM: BOX BOUNDS pp pp pp\n");
    gzputs(fpo, buf);

    sprintf(buf, "%-16.10lf %-16.10lf\n", boxlo[0], boxhi[0]);
    gzputs(fpo, buf);

    sprintf(buf, "%-16.10lf %-16.10lf\n", boxlo[1], boxhi[1]);
    gzputs(fpo, buf);

    sprintf(buf, "%-16.10lf %-16.10lf\n", boxlo[2], boxhi[2]);
    gzputs(fpo, buf);

    sprintf(buf, "ITEM: ATOMS id mol xu yu zu\n");
    gzputs(fpo, buf);

    id = 0;
    for (i = 0; i < nchains; i++) {
        for (j = 0; j < chains[i].n; j++) {
            x = ARRAY2D(chains[i].r, j, 0, 3);
            y = ARRAY2D(chains[i].r, j, 1, 3);
            z = ARRAY2D(chains[i].r, j, 2, 3);
            sprintf(buf, "%d %d %.6f %.6f %.6f\n", id+1, i+1, x, y, z);
            gzputs(fpo, buf);
            id++;
        }
    }

}

int read_frame(gzFile fp, chain_t *chains, int nchains, double *boxlo, double *boxhi)
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
    sscanf(buf, "%lf %lf", &boxlo[0], &boxhi[0]);
    read_line(fp);
    sscanf(buf, "%lf %lf", &boxlo[1], &boxhi[1]);
    read_line(fp);
    sscanf(buf, "%lf %lf", &boxlo[2], &boxhi[2]);

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
        ARRAY2D(chains[mol].r, j, 0, 3) = x;
        ARRAY2D(chains[mol].r, j, 1, 3) = y;
        ARRAY2D(chains[mol].r, j, 2, 3) = z;
        chains[mol].m++;
        assert(chains[mol].m <= chains[mol].n);
    }

    for (i = 0; i < nchains; i++) {
        chains[i].m = 0;
    }

    return step;
}

int read_natoms(char *fn)
{
    int n;
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
        if (strstr(buf, " atoms")) {
            sscanf(buf, "%d", &n);
            fclose(fp);
            return n;
        }
    }

    fclose(fp);
    return -1;
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

double distance(double x1, double y1, double z1, double x2, double y2, double z2)
{
    return sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2));
}

double distance_sq(double x1, double y1, double z1, double x2, double y2, double z2)
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

void print_matrix(double *A, int M, int N)
{
    int i, j;

    for (i = 0; i < M; i++) {
        for (j = 0; j < N; j++) {
            printf("% -6.4f%s", ARRAY2D(A, i, j, N), (j == N-1 ? "\n" : "  "));
        }
    }
}

void transpose(double *A, double *B, int M, int N)
{
    int i, j;

    for (i = 0; i < M; i++) {
        for (j = 0; j < N; j++) {
            ARRAY2D(B, j, i, M) = ARRAY2D(A, i, j, N);
        }
    }
}

void multiply(double *A, double *B, double *C, int M, int N, int K)
{
    int i, j, k;
    double s;

    for (i = 0; i < M; i++) {
        for (j = 0; j < N; j++) {
            s = 0;
            for (k = 0; k < K; k++) {
                s += ARRAY2D(A, i, k, K)*ARRAY2D(B, k, j, N);
            }
            ARRAY2D(C, i, j, N) = s;
        }
    }
}

double det3x3rmo(double *A)
{
    return (A[0]*A[4]*A[8] + A[1]*A[5]*A[6] + A[2]*A[3]*A[7] - A[2]*A[4]*A[6] - A[1]*A[3]*A[8] - A[0]*A[5]*A[7]);
}

void kabsch(chain_t *c, chain_t *ref)
{
    int i, j;
    int N;
    double *P;
    double *Pt;
    double *Q;
    double *A;

    double *SIGMA;
    double *U;
    double *Vt;
    double superb[10];

    double *Ut;
    double *V;
    double d;
    double M[9] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };

    double *R;

    assert(c->n == ref->n);
    N = c->n;

    P  = malloc(N * 3 * sizeof(double));
    Pt = malloc(N * 3 * sizeof(double));
    Q  = malloc(N * 3 * sizeof(double));
    A  = malloc(3 * 3 * sizeof(double));

    assert(P  != NULL);
    assert(Pt != NULL);
    assert(Q  != NULL);
    assert(A  != NULL);

    for (i = 0; i < N; i++) {
        for (j = 0; j < 3; j++) {
            ARRAY2D(P, i, j, 3) = ARRAY2D(c->r, i, j, 3);
            ARRAY2D(Q, i, j, 3) = ARRAY2D(ref->r, i, j, 3);
        }
    }

    LAPACKE_dge_trans(LAPACK_ROW_MAJOR, N, 3, P, 3, Pt, N);
    //transpose(P, Pt, N, 3);

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 3, N, 1.0, Pt, N, Q, 3, 0, A, 3);
    //multiply(Pt, Q, A, 3, 3, N);

    SIGMA = malloc(3 * 3 * sizeof(double));
    U = malloc(3 * 3 * sizeof(double));
    Vt = malloc(3 * 3 * sizeof(double));

    //printf("--A--\n");
    //print_matrix(A, 3, 3);
    //printf("---\n");

    LAPACKE_dgesvd(LAPACK_ROW_MAJOR, 'A', 'A', 3, 3, A, 3, SIGMA, U, 3, Vt, 3, superb);

    //printf("**SIGMA**\n");
    //print_matrix(SIGMA, 3, 3);

    //printf("==U==\n");
    //print_matrix(U, 3, 3);

    //printf("==Vt==\n");
    //print_matrix(Vt, 3, 3);

    Ut = malloc(3 * 3 * sizeof(double));
    V  = malloc(3 * 3 * sizeof(double));

    assert(Ut != NULL);
    assert(V  != NULL);

    LAPACKE_dge_trans(LAPACK_ROW_MAJOR, 3, 3, U, 3, Ut, 3);
    LAPACKE_dge_trans(LAPACK_ROW_MAJOR, 3, 3, Vt, 3, V, 3);

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, 1.0, V, 3, Ut, 3, 0, A, 3);


    d = det3x3rmo(A);

    M[8] = d;
    //printf("---M---\n");
    //print_matrix(M, 3, 3);

    R = malloc(3 * 3 * sizeof(double));
    assert(R != NULL);

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, 1.0, V, 3, M, 3, 0, A, 3);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, 1.0, A, 3, Ut, 3, 0, R, 3);

    //printf("~~~~~~~~~~~~~~~~~~~~~~~~~~ R ~~~~~~~~~~~~~~~~~~~~~~~~~\n");
    //print_matrix(R, 3, 3);
    //printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, N, 3, 1.0, R, 3, Pt, N, 0, Q, N);

    LAPACKE_dge_trans(LAPACK_ROW_MAJOR, 3, N, Q, N, P, 3);

    //printf("________P FINAL_______\n");
    //print_matrix(P, N, 3);
    //printf("^^^^^^^^^^^^^^^^^^^^^^\n");

    /* replace the old coordinates */
    for (i = 0; i < N; i++) {
        for (j = 0; j < 3; j++) {
            ARRAY2D(c->r, i, j, 3) = ARRAY2D(P, i, j, 3);
        }
    }
}

int main(int argc, char *argv[])
{
    int i, j, k;
    int step;
    chain_t chains[NCHAINS];
    chain_t ref[NCHAINS];
    chain_t mean[NCHAINS];
    double rg_min[NCHAINS];
    double clen = CHAINLEN;
    double mass;
    gzFile fp, fpo;
    double xcom, ycom, zcom;
    double rg;
    double xj, yj, zj;
    double xk, yk, zk;
    double boxlo[3], boxhi[3];
    double *sigma, *sigmap, *lambda, *M, *A;
    double omega, Sho;
    double kB = 1.0;
    double T = 2.0;
    double hbar;
    int m;
    int nconf = 0;
    int natoms;

    hbar = 0.18292026/(2.0*M_PI);
    m = 3*clen;
    
    /* consdereing each molecule individually */

    for (i = 0; i < NCHAINS; i++) {
        /* main chains vector */
        chains[i].n = clen;
        chains[i].m = 0;
        chains[i].r = malloc(clen * 3 * sizeof(double));
        assert(chains[i].r != NULL);
        /* reference chains vector */
        ref[i].n = clen;
        ref[i].m = 0;
        ref[i].r = malloc(clen * 3 * sizeof(double));
        assert(ref[i].r != NULL);
        /* mean */
        mean[i].n = clen;
        mean[i].m = 0;
        mean[i].r = malloc(clen * 3 * sizeof(double));
        assert(mean[i].r != NULL);
        for (j = 0; j < clen; j++) {
            ARRAY2D(mean[i].r, j, 0, 3) = 0;
            ARRAY2D(mean[i].r, j, 1, 3) = 0;
            ARRAY2D(mean[i].r, j, 2, 3) = 0;
        }
        /* min RG found for each chain */
        rg_min[i] = DBL_MAX;
    }

    mass = read_mass("data.lmp");
    assert(mass > 0);

    natoms = read_natoms("data.lmp");
    assert(natoms > 0);

#if 0
    fp = gzopen("traj.dump.gz", "r");
    assert(fp != NULL);

    if (gzbuffer(fp, 0xF0000) == -1) {
        fprintf(stderr, "error: gzbuffer\n");
        return EXIT_FAILURE;
    }

    fpo = gzopen("temp.dump.gz", "w");
    assert(fpo != NULL);

    if (gzbuffer(fpo, 0xF0000) == -1) {
        fprintf(stderr, "error: gzbuffer\n");
        return EXIT_FAILURE;
    }



    while ((step = read_frame(fp, chains, NCHAINS, boxlo, boxhi)) >= 0) {

        printf("1ST PASS: %d\n", step);

        /* Remove COM from each molecule */
        for (i = 0; i < NCHAINS; i++) {
            xcom = 0;
            ycom = 0;
            zcom = 0;
            for (j = 0; j < chains[i].n; j++) {
                xcom += (ARRAY2D(chains[i].r, j, 0, 3) * mass);
                ycom += (ARRAY2D(chains[i].r, j, 1, 3) * mass);
                zcom += (ARRAY2D(chains[i].r, j, 2, 3) * mass);
            }
            xcom /= (chains[i].n * mass);
            ycom /= (chains[i].n * mass);
            zcom /= (chains[i].n * mass);
            for (j = 0; j < chains[i].n; j++) {
                ARRAY2D(chains[i].r, j, 0, 3) -= xcom;
                ARRAY2D(chains[i].r, j, 1, 3) -= ycom;
                ARRAY2D(chains[i].r, j, 2, 3) -= zcom;
            }
        }

        /* Save reference */
        for (i = 0; i < NCHAINS; i++) {
            rg = 0;
            for (j = 0; j < chains[i].n; j++) {
                xj = ARRAY2D(chains[i].r, j, 0, 3);
                yj = ARRAY2D(chains[i].r, j, 1, 3);
                zj = ARRAY2D(chains[i].r, j, 2, 3);
                for (k = 0; k < chains[i].n; k++) {
                    xk = ARRAY2D(chains[j].r, k, 0, 3);
                    yk = ARRAY2D(chains[i].r, k, 1, 3);
                    zk = ARRAY2D(chains[i].r, k, 2, 3);
                    rg += distance_sq(xj, yj, zj, xk, yk, zk);
                }
                rg /= (2.0 * chains[i].n * chains[i].n);
                rg = sqrt(rg);
            }
            if (rg < rg_min[i]) {
                rg_min[i] = rg;
                copy_chain(&chains[i], &ref[i]);
            }
        }

        /* Remove rigid-body rotation */
        for (i = 0; i < NCHAINS; i++) {
            kabsch(&chains[i], &ref[i]);
        }

        nconf++;

        /* Write frame */
        write_frame(fpo, chains, step, NCHAINS, natoms, boxlo, boxhi);

        /* Compute fluctuations */
        for (i = 0; i < NCHAINS; i++) {
            for (j = 0; j < mean[i].n; j++) {
                ARRAY2D(mean[i].r, j, 0, 3) += ((ARRAY2D(chains[i].r, j, 0, 3) - ARRAY2D(mean[i].r, j, 0, 3))/((double)nconf));
                ARRAY2D(mean[i].r, j, 1, 3) += ((ARRAY2D(chains[i].r, j, 1, 3) - ARRAY2D(mean[i].r, j, 1, 3))/((double)nconf));
                ARRAY2D(mean[i].r, j, 2, 3) += ((ARRAY2D(chains[i].r, j, 2, 3) - ARRAY2D(mean[i].r, j, 2, 3))/((double)nconf));
            }
        }

        //break;
    }

    gzclose(fp);
    gzclose(fpo);
#endif 
    
    sigma = malloc(m * m * sizeof(double));
    sigmap = malloc(m * m * sizeof(double));
    lambda = malloc(m * sizeof(double));
    M = malloc(m * m * sizeof(double));
    A = malloc(m * m * sizeof(double));
    assert(sigma != NULL);
    assert(sigmap != NULL);
    assert(M != NULL);
    assert(A != NULL);
    assert(lambda != NULL);
    m = 3*clen;

    for (i = 0; i < m; i++) {
        for (j = 0; j < m; j++) {
            ARRAY2D(sigma, i, j, m) = 0;
        }
    }

    fp = gzopen("temp.dump.gz", "r");
    assert(fp != NULL);

    if (gzbuffer(fp, 0xF0000) == -1) {
        fprintf(stderr, "error: gzbuffer\n");
        return EXIT_FAILURE;
    }

    nconf = 0;
    while ((step = read_frame(fp, chains, NCHAINS, boxlo, boxhi)) >= 0) {
        printf("2ND PASS: %d\n", step);
        nconf++;
        for (i = 0; i < NCHAINS; i++) {
            for (j = 0; j < m; j++) {
                for (k = 0; k < m; k++) {
                    double sigma_new = (chains[i].r[j] - mean[i].r[j])*(chains[i].r[k] - mean[i].r[k]);
                    ARRAY2D(sigma, j, k, m) += ((sigma_new - ARRAY2D(sigma, j, k, m)))/((double)nconf);
                }
            }
        }
    }

    fprintf(stderr, "=> CONFIGURATIONS PARSED: %d\n", nconf);


#if 1
    for (i = 0; i < m; i++) {
        for (j = 0; j < m; j++) {
            if (i == j) {
                ARRAY2D(M, i, j, m) = sqrt(m);
            } else {
                ARRAY2D(M, i, j, m) = 0;
            }
        }
    }
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, m, m, 1.0, M, m, sigma, m, 0, A, m);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, m, m, 1.0, A, m, M, m, 0, sigmap, m);
    for (i = 0; i < m; i++) {
        for (j = 0; j < m; j++) {
            if (i > j) {
                ARRAY2D(sigmap, i, j, m) = 0;
            }
        }
    }
    LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'N', 'U', m, sigmap, m, lambda);
#endif

    Sho = 0;
    for (i = 0; i < 6; i++) {
        fprintf(stderr, "%20.10f  ", lambda[i]);
    }
    fprintf(stderr, "\n");
            
    return 0;
}
