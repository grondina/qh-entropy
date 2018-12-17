#include <assert.h>
#include <cblas.h>
#include <lapacke.h>
#include <lapacke_utils.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>

#include "data.h"
#include "io.h"
#include "traj.h"
#include "util.h"

void init_frame(struct frame *frame, struct data *data)
{
    if (frame == NULL)
        return;

    frame->nmols = data->nmols;
    frame->mol = malloc(frame->nmols * (sizeof (struct molecule)));
    assert(frame->mol != NULL);

    for (int i = 0; i < frame->nmols; ++i) {
        init_molecule(&frame->mol[i], data->molsize);
    }

    frame->step = -1;
}

void free_frame(struct frame *frame)
{
    if (frame == NULL)
        return;

    if (frame->mol == NULL)
        return;

    for (int i = 0; i < frame->nmols; ++i)
        free_molecule(&frame->mol[i]);

    free(frame->mol);
}

void parse_pass1(const char *fndump, struct data *data, struct molecule *refmols)
{
    struct frame frame;
    init_frame(&frame, data);

    assert(fndump != NULL);
    gzFile fpdump = gzopen(fndump, "r");
    assert(fpdump != NULL);

    if (gzbuffer(fpdump, 0xF0000) == -1) {
        fprintf(stderr, "error: gzbuffer\n");
            exit(EXIT_FAILURE);
    }

    long int step;
    while ((step = read_frame(fpdump, &frame, data)) >= 0) {

        printf("(1) TIMESTEP: %ld\n", step);

        /* Obtaining reference based on RoG (smallest) */
        for (int i = 0; i < frame.nmols; ++i) {
            frame.mol[i].gyr = gyration(&frame.mol[i]);
            if (frame.mol[i].gyr < refmols[i].gyr)
                copy_molecule(&refmols[i], &frame.mol[i]);
        }
    }

    /* Make sure we stopped reading because file ended */
    if (!gzeof(fpdump)) {
        int err;
        fprintf(stderr, "error: %s\n", gzerror(fpdump, &err));
        exit(EXIT_FAILURE);
    }


    /* Remove COM of reference */
    for (int i = 0; i < frame.nmols; ++i)
        remove_com(&refmols[i]);

    gzclose(fpdump);
    free_frame(&frame);
}

void parse_pass2(const char *fndump, const char *fntemp, struct data *data, struct molecule *refmols, struct molecule *avemols)
{
    struct frame frame;
    init_frame(&frame, data);

    assert(fndump != NULL);
    gzFile fpdump = gzopen(fndump, "r");
    assert(fpdump != NULL);

    if (gzbuffer(fpdump, 0xF0000) == -1) {
        fprintf(stderr, "error: gzbuffer\n");
        exit(EXIT_FAILURE);
    }

    gzFile fptemp = gzopen(fntemp, "w");
    assert(fptemp != NULL);

    if (gzbuffer(fptemp, 0xF0000) == -1) {
        fprintf(stderr, "error: gzbuffer\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < data->nmols; ++i) {
        for (int j = 0; j < avemols[i].n; ++j) {
            avemols[i].R[j][0] = 0;
            avemols[i].R[j][1] = 0;
            avemols[i].R[j][2] = 0;
        }
    }

    long int step;
    int nframes = 0;

    while ((step = read_frame(fpdump, &frame, data)) >= 0) {

        printf("(2) TIMESTEP: %ld\n", step);
        nframes++;

        for (int i = 0; i < frame.nmols; ++i) {

            /* Remove center of mass */
            remove_com(&frame.mol[i]);

            /* Remove right body rotation */
            kabsch(&frame.mol[i], &refmols[i]);

            /* Update mean positions (cummulative rolling average) */
            for (int j = 0; j < avemols[i].n; ++j) {
                avemols[i].R[j][0] += (frame.mol[i].R[j][0] - avemols[i].R[j][0])/((double)nframes);
                avemols[i].R[j][1] += (frame.mol[i].R[j][1] - avemols[i].R[j][1])/((double)nframes);
                avemols[i].R[j][2] += (frame.mol[i].R[j][2] - avemols[i].R[j][2])/((double)nframes);
            }
        }

        /* Write frame with processed coordinates */
        write_frame_bin(fptemp, &frame, data);
    }

    /* Make sure we stopped reading because file ended */
    if (!gzeof(fpdump)) {
        int err;
        fprintf(stderr, "error: %s\n", gzerror(fpdump, &err));
        exit(EXIT_FAILURE);
    }

    gzclose(fpdump);
    gzclose(fptemp);
    free_frame(&frame);
}

void parse_pass3(const char *fndump, struct data *data, struct molecule *avemols, int n, double (*sigma)[n][n], double (*M)[n][n])
{
    struct frame frame;
    init_frame(&frame, data);

    assert(fndump != NULL);
    gzFile fpdump = gzopen(fndump, "r");
    assert(fpdump != NULL);

    if (gzbuffer(fpdump, 0xF0000) == -1) {
        fprintf(stderr, "error: gzbuffer\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < data->nmols; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                sigma[i][j][k] = 0;
                M[i][j][k] = 0;
            }
        }
    }

    long int step;
    int nframes = 0;

    while ((step = read_frame_bin(fpdump, &frame, data)) >= 0) {

        printf("(3) TIMESTEP: %ld\n", step);
        nframes++;

        if (nframes == 1) {
            /* Build diagonal mass matrix */
            for (int i = 0; i < frame.nmols; ++i) {
                for (int j = 0; j < frame.mol[i].m; ++j) {
                    int k = 3 * j;
                    M[i][k + 0][k + 0] = sqrt(frame.mol[i].mass[j]);
                    M[i][k + 1][k + 1] = sqrt(frame.mol[i].mass[j]);
                    M[i][k + 2][k + 2] = sqrt(frame.mol[i].mass[j]);
                }
            }
        }

        for (int i = 0; i < frame.nmols; ++i) {
            double *x = &(frame.mol[i].R[0][0]);
            double *y = &(  avemols[i].R[0][0]);
            for (int j = 0; j < n; ++j) {
                for (int k = 0; k < n; ++k) {
                    double value = (x[j] - y[j]) * (x[k] - y[k]);
                    sigma[i][j][k] += (value - sigma[i][j][k])/((double)nframes);
                }
            }
        }
    }

    if (!gzeof(fpdump)) {
        int err;
        fprintf(stderr, "error: %s\n", gzerror(fpdump, &err));
        exit(EXIT_FAILURE);
    }

    gzclose(fpdump);
    free_frame(&frame);
}

double *entropy(int n, struct data *data, struct arguments *arguments, double (*sigma)[n][n], double (*M)[n][n])
{
    double kB = arguments->kB;
    double hbar = arguments->h/(2.0 * M_PI);
    double temp = arguments->temp;

    double (*sigmap)[n][n] = malloc(sizeof(double[data->nmols][n][n]));;
    assert(sigmap != NULL);

    double (*X)[n] = malloc(sizeof(double[n][n]));
    assert(X != NULL);

    double *lambda = malloc(sizeof(double[n]));
    assert(lambda != NULL);

    double *S = malloc(data->nmols * sizeof(double));
    assert(S != NULL);

    double *pX = &(X[0][0]);

    for (int i = 0; i < data->nmols; ++i) {

        double *pS  = &( sigma[i][0][0]);
        double *pSp = &(sigmap[i][0][0]);
        double *pM  = &(     M[i][0][0]);

        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, pM, n, pS, n, 0,  pX, n);
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, pX, n, pM, n, 0, pSp, n);

        LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'N', 'U', n, pSp, n, lambda);

        S[i] = 0;
        for (int j = 6; j < n; ++j) {
            assert(fabs(lambda[j]) > 1.0e-10);
            double omega = sqrt((kB * temp)/lambda[j]);
            double factor = (hbar * omega)/(kB * temp);
            S[i] += (factor/(exp(factor) - 1.0) - log(1.0 - exp(-factor)));
            S[i] *= kB;
        }
    }

    free(sigmap);
    free(X);
    free(lambda);

    return S;
}
