#ifndef _DATA_H
#define _DATA_H

struct data {
    int ntypes;
    int natoms;
    int nmols;
    int molsize;
    int *type;
    double temp;
    double xlo, xhi;
    double ylo, yhi;
    double zlo, zhi;
    double xlen;        /* box length in x direction */
    double ylen;        /* box length in y direction */
    double zlen;        /* box length in z direction */
    double *mass;
};

struct molecule {
    int id;         /* LAMMPS id of the molecule */
    int n;          /* max number of atoms in this molecule */
    int m;          /* number of atoms already read */
    int *atoms;     /* sorted array of n atom indices   */
    int *types;     /* types of atoms */
    double gyr;     /* radius of gyration */
    double *mass;
    double (*R)[3]; /* coordinates, n rows x 3 columns  */
};

extern void init_data(struct data *data);
extern void free_data(struct data *data);

extern void init_molecule(struct molecule *molecule, int n);
extern void free_molecule(struct molecule *molecule);

extern void print_atoms(struct molecule *mol);

extern double get_mass(struct data *data, int type);

extern void copy_molecule(struct molecule *m2, struct molecule *m1);

extern struct molecule *init_molecule_array(struct data *data);
extern void free_molecule_array(struct molecule *array, struct data *data);



#endif
