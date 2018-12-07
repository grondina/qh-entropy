#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "io.h"
#include "data.h"

static char buf[BUFSIZ];

int read_data(char *fndata, struct data *data)
{
    int type;
    int ntypes;
    double mass;
    FILE *fp;

    assert(fndata != NULL);
    assert(data != NULL);

    fp = fopen(fndata, "r");
    if (fp == NULL) {
        perror("read_data");
        return -1;
    }

    while (fgets(buf, BUFSIZ, fp)) {

        if (strstr(buf, "atom types")) {

            sscanf(buf, "%d", &ntypes);
            data->ntypes = ntypes;

            continue;
        }

        if (strstr(buf, "Masses")) {

            /* Make sure we know how many types are there */
            assert(data->ntypes > 0);

            /* Allocate vectors */
            data->mass = malloc(data->ntypes * sizeof(double));
            data->type = malloc(data->ntypes * sizeof(double));
            assert(data->mass != NULL);
            assert(data->type != NULL);

            /* Skip blank line */
            assert(fgets(buf, BUFSIZ, fp));

            /* Read ntypes masses */
            for (int i = 0; i < data->ntypes; ++i) {

                assert(fgets(buf, BUFSIZ, fp));

                if (sscanf(buf, "%d %lf", &type, &mass) == 2) {
                    assert(type > 0);
                    data->type[i] = type;
                    data->mass[i] = mass;
                } else {
                    fprintf(stderr, "read_data: not enough masses\n");
                    return -1;
                }
            }

            break;
        }
    }

    if (ferror(fp)) {
        perror("read_data");
        return -1;
    }

    fclose(fp);
    return 0;
}

int read_frame(void)
{
    return 0;
}

void write_frame(void)
{
}
