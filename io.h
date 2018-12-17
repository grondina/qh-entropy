#ifndef _IO_H
#define _IO_H

#include "data.h"
#include "traj.h"

extern int read_data(char *fndata, struct data *data);
extern long int read_frame(gzFile fp, struct frame *frame, struct data *data);
extern void write_frame(gzFile fp, struct frame *frame, struct data *data);
extern void write_entropy(char *fn, struct data *data, double *S);
extern void init_buf(struct data *data);
extern void free_buf(void);

void write_frame_bin(gzFile fp, struct frame *frame, struct data *data);
long int read_frame_bin(gzFile fp, struct frame *frame, struct data *data);

#endif
