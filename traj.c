#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include "data.h"
#include "io.h"
#include "ref.h"
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

void parse_pass2(const char *fndump, const char *fntemp, struct data *data, struct molecule *refmols)
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

    long int step;
    while ((step = read_frame(fpdump, &frame, data)) >= 0) {

        printf("(2) TIMESTEP: %ld\n", step);

        for (int i = 0; i < frame.nmols; ++i) {

            /* Remove center of mass */
            remove_com(&frame.mol[i]);

            /* Remove right body rotation */
            kabsch(&frame.mol[i], &refmols[i]);
        }

        /* Write frame with processed coordinates */
        write_frame(fptemp, &frame, data);
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
