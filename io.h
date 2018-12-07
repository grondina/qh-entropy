#ifndef _IO_H
#define _IO_H

#include "data.h"

extern int read_data(char *fndata, struct data *data);
extern int read_frame(void);
extern void write_frame(void);

#endif
