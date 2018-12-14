#include <assert.h>
#include <cblas.h>
#include <float.h>
#include <lapacke.h>
#include <lapacke_utils.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>

#include "data.h"
#include "io.h"
#include "par.h"
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

void parse_pass1(const char *fndump, struct data *data, struct molecule *refmols, int nthreads)
{
    assert(fndump != NULL);
    gzFile fpdump = gzopen(fndump, "r");
    assert(fpdump != NULL);

    if (gzbuffer(fpdump, 0xF0000) == -1) {
        fprintf(stderr, "error: gzbuffer\n");
            exit(EXIT_FAILURE);
    }

    long int *steps = malloc(nthreads * sizeof(long int));
    assert(steps != NULL);

    struct thread *threads = init_threads(nthreads, data);

    struct molecule **refmols_t = malloc(nthreads * sizeof(struct molecule *));
    assert(refmols_t != NULL);
    for (int i = 0; i < nthreads; ++i)
        refmols_t[i] = init_molecule_array(data);

    while (read_frames(fpdump, data, threads, nthreads, steps) > 0) {
        printf("(1) TIMESTEPS: ");
        for (int t = 0; t < nthreads; ++t)
            if (steps[t] >= 0) printf("%10ld ", steps[t]);
        printf("\n");
        #pragma omp parallel default(none) shared(steps, threads, refmols_t, stderr) num_threads(nthreads)
        {
            int tid = omp_get_thread_num();
            if (steps[tid] >= 0) {
                //fprintf(stderr, "hello from thread %d\n", tid);
                struct frame *frame = threads[tid].frame;
                for (int i = 0; i < frame->nmols; ++i) {
                    frame->mol[i].gyr = gyration(&frame->mol[i]);
                    if (frame->mol[i].gyr < refmols_t[tid][i].gyr)
                        copy_molecule(&refmols_t[tid][i], &frame->mol[i]);
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

    /* Take the reference with lowest RoG */
    for (int i = 0; i < data->nmols; ++i) {
        int idxlow = -1;
        double gyrlow = DBL_MAX;
        for (int t = 0; t < nthreads; ++t) {
            if (refmols_t[t][i].gyr < gyrlow) {
                idxlow = t;
                gyrlow = refmols_t[t][i].gyr;
            }
        }
        assert(idxlow >= 0);
        copy_molecule(&refmols[i], &refmols_t[idxlow][i]);
    }

    /* Remove COM of reference */
    for (int i = 0; i < data->nmols; ++i)
        remove_com(&refmols[i]);

    free(steps);
    free_threads(threads, nthreads);
    for (int i = 0; i < nthreads; ++i)
        free_molecule_array(refmols_t[i], data);
    free(refmols_t);
}

void parse_pass2(const char *fndump, const char *fntemp, struct data *data, struct molecule *refmols, struct molecule *avemols, int nthreads)
{
    assert(fndump != NULL);
    gzFile fpdump = gzopen(fndump, "r");
    assert(fpdump != NULL);

    if (gzbuffer(fpdump, 0xF0000) == -1) {
        fprintf(stderr, "error: gzbuffer\n");
            exit(EXIT_FAILURE);
    }

    assert(fntemp != NULL);
    gzFile fptemp = gzopen(fntemp, "w");
    assert(fptemp != NULL);

    if (gzbuffer(fptemp, 0xF0000) == -1) {
        fprintf(stderr, "error: gzbuffer\n");
        exit(EXIT_FAILURE);
    }


    long int *steps = malloc(nthreads * sizeof(long int));
    assert(steps != NULL);

    struct thread *threads = init_threads(nthreads, data);

    struct molecule **avemols_t = malloc(nthreads * sizeof(struct molecule *));
    assert(avemols_t != NULL);
    for (int i = 0; i < nthreads; ++i)
        avemols_t[i] = init_molecule_array(data);

    for (int t = 0; t < nthreads; ++t) {
        for (int i = 0; i < data->nmols; ++i) {
            for (int j = 0; j < avemols[i].n; ++j) {
                avemols_t[t][i].R[j][0] = 0;
                avemols_t[t][i].R[j][1] = 0;
                avemols_t[t][i].R[j][2] = 0;
            }
        }
    }

    int framecnt;
    int nframes = 0;

    while ((framecnt = read_frames(fpdump, data, threads, nthreads, steps)) > 0) {
        nframes += framecnt;
        printf("(2) TIMESTEPS: ");
        for (int t = 0; t < nthreads; ++t)
            if (steps[t] >= 0) printf("%10ld ", steps[t]);
        printf("\n");
        #pragma omp parallel default(none) shared(steps, threads, avemols_t, refmols, stderr) num_threads(nthreads)
        {
            int tid = omp_get_thread_num();
            if (steps[tid] >= 0) {
                struct frame *frame = threads[tid].frame;
                for (int i = 0; i < frame->nmols; ++i) {
                    remove_com(&frame->mol[i]);
                    kabsch(&frame->mol[i], &refmols[i]);
                    for (int j = 0; j < avemols_t[tid][i].n; ++j) {
                        avemols_t[tid][i].R[j][0] += frame->mol[i].R[j][0];
                        avemols_t[tid][i].R[j][1] += frame->mol[i].R[j][1];
                        avemols_t[tid][i].R[j][2] += frame->mol[i].R[j][2];
                    }
                }
            }
        }
        for (int t = 0; t < nthreads; ++t)
            if (steps[t] >= 0) 
                write_frame(fptemp, threads[t].frame, data);
    }

    if (!gzeof(fpdump)) {
        int err;
        fprintf(stderr, "error: %s\n", gzerror(fpdump, &err));
        exit(EXIT_FAILURE);
    }
    gzclose(fpdump);
    gzclose(fptemp);

    /* Reduce sum of coordinates and obtain average */
    for (int i = 0; i < data->nmols; ++i) {
        for (int j = 0; j < avemols[i].n; ++j) {
            avemols[i].R[j][0] = 0;
            avemols[i].R[j][1] = 0;
            avemols[i].R[j][2] = 0;
        }
    }
    for (int i = 0; i < data->nmols; ++i) {
        for (int t = 0; t < nthreads; ++t) {
            for (int j = 0; j < avemols[i].n; ++j) {
                avemols[i].R[j][0] += avemols_t[t][i].R[j][0];
                avemols[i].R[j][1] += avemols_t[t][i].R[j][1];
                avemols[i].R[j][2] += avemols_t[t][i].R[j][2];
            }
        }
    }
    for (int i = 0; i < data->nmols; ++i) {
        for (int j = 0; j < avemols[i].n; ++j) {
            avemols[i].R[j][0] /= (double)nframes;
            avemols[i].R[j][1] /= (double)nframes;
            avemols[i].R[j][2] /= (double)nframes;
        }
    }

    free(steps);
    free_threads(threads, nthreads);
    for (int i = 0; i < nthreads; ++i)
        free_molecule_array(avemols_t[i], data);
    free(avemols_t);
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

    while ((step = read_frame(fpdump, &frame, data)) >= 0) {

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
