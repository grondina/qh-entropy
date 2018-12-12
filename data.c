#include <assert.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "data.h"

void copy_molecule(struct molecule *m2, struct molecule *m1)
{
    if (m1 == NULL)
        return;

    if (m2 == NULL)
        return;

    m2->n = m1->n;
    m2->m = m1->m;
    m2->gyr = m1->gyr;

    assert(m1->atoms != NULL);
    assert(m2->atoms != NULL);
    memcpy(m2->atoms, m1->atoms, (m1->n * sizeof(int)));

    assert(m1->types != NULL);
    assert(m2->types != NULL);
    memcpy(m2->types, m1->types, (m1->n * sizeof(int)));

    assert(m1->mass != NULL);
    assert(m2->mass != NULL);
    memcpy(m2->mass, m1->mass, (m1->n * sizeof(double)));

    assert(m1->R != NULL);
    assert(m2->R != NULL);
    memcpy(m2->R, m1->R, (m1->n * 3 * sizeof(double)));
}

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
    data->xlo = -1;
    data->xhi = -1;
    data->ylo = -1;
    data->yhi = -1;
    data->zlo = -1;
    data->zhi = -1;
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
    mol->id = -1;
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

struct molecule *init_molecule_array(struct data *data)
{
    struct molecule *array = malloc(data->nmols * sizeof(struct molecule));
    assert(array != NULL);

    for (int i = 0; i < data->nmols; ++i)
        init_molecule(&array[i], data->molsize);

    return array;
}

void free_molecule_array(struct molecule *array, struct data *data)
{
    if (array == NULL)
        return;

    if (data == NULL)
        return;

    for (int i = 0; i < data->nmols; ++i)
        free_molecule(&array[i]);

    free(array)
}
