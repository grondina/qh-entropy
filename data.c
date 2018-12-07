#include <stdio.h>
#include <stdlib.h>
#include "data.h"

void init_data(struct data *data)
{
    data->ntypes = 0;
    data->mass = NULL;
    data->type = NULL;
}

void free_data(struct data *data)
{
    free(data->mass);
    free(data->type);
}
