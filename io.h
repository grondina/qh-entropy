#ifndef _IO_H
#define _IO_H

#include "data.h"
#include "traj.h"

extern int read_data(char *fndata, struct data *data);
extern long int read_frame(gzFile fp, struct frame *frame, struct data *data);
extern void write_frame(gzFile fp, struct frame *frame, struct data *data);
#endif
