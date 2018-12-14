#ifndef _IO_H
#define _IO_H

#include "data.h"
#include "par.h"
#include "traj.h"

extern int read_data(char *fndata, struct data *data);
extern int read_frames(gzFile fp, struct data *data, struct thread *threads, int nframes, long int *ret);
extern long int read_frame(gzFile fp, struct frame *frame, struct data *data);
extern void write_frame(gzFile fp, struct frame *frame, struct data *data);
extern void write_entropy(char *fn, struct data *data, double *S);

#endif
