#ifndef _DATA_H
#define _DATA_H

struct data {
    int ntypes;
    int natoms;
    int nmol;
    int *type;
    double xlen;        /* box length in x direction */
    double ylen;        /* box length in y direction */
    double zlen;        /* box length in z direction */
    double *mass;
};

extern void init_data(struct data *data);
extern void free_data(struct data *data);

#endif
