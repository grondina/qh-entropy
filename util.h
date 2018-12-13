#ifndef _UTIL_H
#define _UTIL_H

#include "data.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

struct arguments {
    int n;              /* chain length */
    double kB;          /* Boltzmann constant */
    double h;           /* Planck's constant */
    double temp;        /* temperature */
    char *fndata;       /* name of data file */
    char *fndump;       /* name of dump file */
    char *fntemp;       /* name of temp dump file */
};

extern int cmpint(const void *p1, const void *p2);
extern double gyration(struct molecule *mol);
extern void kabsch(struct molecule *mol, struct molecule *ref);
extern void remove_com(struct molecule *mol);
extern void print_matrix(int N, int M, double (*A)[M]);
extern void print_vector(int N, double *v);

#endif
