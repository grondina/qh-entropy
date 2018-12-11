#ifndef _DATA_H
#define _DATA_H

struct data {
    int ntypes;
    int natoms;
    int nmols;
    int molsize;
    int *type;
    double temp;
    double xlen;        /* box length in x direction */
    double ylen;        /* box length in y direction */
    double zlen;        /* box length in z direction */
    double *mass;
};

struct molecule {
    int n;          /* max number of atoms in this molecule */
    int m;          /* number of atoms already read */
    int *atoms;     /* sorted array of n atom indices   */
    int *types;     /* types of atoms */
    double *mass;
    double (*R)[3]; /* coordinates, n rows x 3 columns  */
};

extern void init_data(struct data *data);
extern void free_data(struct data *data);

extern void init_molecule(struct molecule *molecule, int n);
extern void free_molecule(struct molecule *molecule);

extern void print_atoms(struct molecule *mol);

double get_mass(struct data *data, int type);

#endif
