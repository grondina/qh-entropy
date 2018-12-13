#ifndef _TRAJ_H
#define _TRAJ_H

#include "data.h"
#include "util.h"

struct frame {
    long int step;
    int nmols;
    struct molecule *mol;
};

extern void init_frame(struct frame *frame, struct data *data);
extern void free_frame(struct frame *frame);
extern void parse_pass1(const char *fndump, struct data *data, struct molecule *refmols);
extern void parse_pass2(const char *fndump, const char *fntemp, struct data *data, struct molecule *refmols, struct molecule *avemols);
extern void parse_pass3(const char *fndump, struct data *data, struct molecule *avemols, int n, double (*sigma)[n][n], double (*M)[n][n]);
extern double *entropy(int n, struct data *data, struct arguments *arguments, double (*sigma)[n][n], double (*M)[n][n]);

#endif
