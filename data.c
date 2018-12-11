#include <assert.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include "data.h"

void init_data(struct data *data)
{
    if (data == NULL)
        return;

    data->ntypes = -1;
    data->ntypes = -1;
    data->nmols= -1;
    data->molsize = -1;
    data->temp = -1;
    data->mass = NULL;
    data->type = NULL;
    data->xlen = -1;
    data->ylen = -1;
    data->zlen = -1;
}

void free_data(struct data *data)
{
    if (data == NULL)
        return;

    free(data->mass);
    free(data->type);
}

void init_molecule(struct molecule *mol, int n)
{
    if (mol == NULL)
        return;

    mol->m = 0;
    mol->n = n;
    mol->gyr = DBL_MAX;

    mol->atoms = malloc(mol->n * (sizeof (int)));
    assert(mol->atoms != NULL);

    mol->types = malloc(mol->n * (sizeof(int)));
    assert(mol->types != NULL);

    mol->mass = malloc(mol->n * (sizeof(double)));
    assert(mol->mass != NULL);

    mol->R = malloc(sizeof (double([mol->n][3])));
    assert(mol->R != NULL);
}


void free_molecule(struct molecule *molecule)
{
    if (molecule == NULL)
        return;

    free(molecule->atoms);
    free(molecule->types);
    free(molecule->mass);
    free(molecule->R);
}

void print_atoms(struct molecule *mol)
{
    if (mol == NULL)
        return;

    for (int i = 0; i < mol->m; ++i) {
        char *s = (i == (mol->m - 1)) ? "\n" : "   ";
        printf("%3d%s", mol->atoms[i], s);
    }

    printf("-----------------\n");
    for (int i = 0; i < mol->m; ++i) {
        printf("%10.6f    %10.6f    %10.6f   (%3d, %4.2f)\n",
              mol->R[i][0], mol->R[i][1], mol->R[i][2],
              mol->types[i], mol->mass[i]);
    }
}

double get_mass(struct data *data, int type)
{
    if (data == NULL)
        return -1;

    for (int i = 0; i < data->ntypes; ++i) {
        if (type == data->type[i])
            return data->mass[i];
    }

    return -1;
}
