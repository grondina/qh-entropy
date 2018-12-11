#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "data.h"
#include "ref.h"
#include "traj.h"

struct molecule *init_reference(struct data *data)
{
    struct molecule *refmols = malloc(data->nmols * sizeof(struct molecule));
    assert(refmols != NULL);

    for (int i = 0; i < data->nmols; ++i)
        init_molecule(&refmols[i], data->molsize);

    return refmols;
}

void free_reference(struct molecule *refmols, struct data *data)
{
    if (refmols == NULL)
        return;

    for (int i = 0; i < data->nmols; ++i)
        free_molecule(&refmols[i]);

    free(refmols);
}

void print_reference(struct molecule *refmols, struct data *data)
{
    for (int i = 0; i < data->nmols; ++i)
        printf("mol %3d  :  RoG = %6.2f\n", i, refmols[i].gyr);
}
