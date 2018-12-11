#include <assert.h>
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
    mol->atoms = malloc(mol->n * (sizeof (int)));
    assert(mol->atoms != NULL);

    mol->R = malloc(sizeof (double([mol->n][3])));
    assert(mol->R != NULL);
}


void free_molecule(struct molecule *molecule)
{
    if (molecule == NULL)
        return;

    free(molecule->atoms);
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
        printf("%10.6f    %10.6f    %10.6f\n",
              mol->R[i][0], mol->R[i][1], mol->R[i][2]);
    }
}
