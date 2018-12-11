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
