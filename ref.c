#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "data.h"
#include "ref.h"
#include "traj.h"

struct molecule *refmols = NULL;

void init_reference(struct data *data)
{
    refmols = malloc(data->nmols * sizeof(struct molecule));
    assert(refmols != NULL);

    for (int i = 0; i < data->nmols; ++i)
        init_molecule(&refmols[i], data->molsize);
}

void free_reference(struct data *data)
{
    if (refmols == NULL)
        return;

    for (int i = 0; i < data->nmols; ++i)
        free_molecule(&refmols[i]);

    free(refmols);
}
