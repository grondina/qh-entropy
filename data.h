#ifndef _DATA_H
#define _DATA_H

struct data {
    int ntypes;
    int *type;
    double *mass;
};

extern void init_data(struct data *data);
extern void free_data(struct data *data);

#endif
