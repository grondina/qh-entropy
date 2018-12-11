#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include "data.h"
#include "io.h"
#include "traj.h"

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

void parse_pass1(const char *fndump, struct data *data)
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

        printf("TIMESTEP: %ld\n", step);

#if 0
        for (int i = 0; i < frame.nmols; ++i) {
            printf("molecule %3d:   ", i);
            print_atoms(&frame.mol[i]);
        }
#endif

        double xcom;
        double ycom;
        double zcom;
        double mass;

        for (int i = 0; i < frame.nmols; ++i) {
            xcom = 0;
            ycom = 0;
            zcom = 0;
            mass = 0;
            for (int j = 0; j < frame.mol[i].m; ++j) {
                xcom += (frame.mol[i].R[j][0] * frame.mol[i].mass[j]);
                ycom += (frame.mol[i].R[j][1] * frame.mol[i].mass[j]);
                zcom += (frame.mol[i].R[j][2] * frame.mol[i].mass[j]);
                mass += frame.mol[i].mass[j];
            }
            xcom /= mass;
            ycom /= mass;
            zcom /= mass;
            for (int j = 0; j < frame.mol[i].m; ++j) {
                frame.mol[i].R[j][0] -= xcom;
                frame.mol[i].R[j][1] -= xcom;
                frame.mol[i].R[j][2] -= xcom;
            }
        }
    }

    gzclose(fpdump);
    free_frame(&frame);
}
