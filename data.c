#include <stdio.h>
#include <stdlib.h>
#include "data.h"

void init_data(struct data *data)
{
    data->ntypes = -1;
    data->ntypes = -1;
    data->nmol= -1;
    data->mass = NULL;
    data->type = NULL;
    data->xlen = -1;
    data->ylen = -1;
    data->zlen = -1;
}

void free_data(struct data *data)
{
    free(data->mass);
    free(data->type);
}
