#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "data.h"
#include "par.h"
#include "traj.h"

struct thread *init_threads(int nthreads, struct data *data)
{
    struct thread *threads = malloc(nthreads * sizeof(struct thread));
    assert(threads != NULL);

    for (int i = 0; i < nthreads; ++i) {
        threads[i].frame = malloc(sizeof(struct frame));
        assert(threads[i].frame != NULL);
        init_frame(threads[i].frame, data);
    }

    return threads;
}

void free_threads(struct thread *threads, int nthreads)
{
    if (threads == NULL)
        return;

    for (int i = 0; i < nthreads; ++i)
        free(threads[i].frame);

    free(threads);
}
