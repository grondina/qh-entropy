#ifndef _PAR_H
#define _PAR_H

#include "data.h"
#include "traj.h"

struct thread {
    struct frame *frame;
};

extern struct thread *init_threads(int nthreads, struct data *data);
extern void free_threads(struct thread *threads, int nthreads);

#endif
