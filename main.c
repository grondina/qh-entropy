#include <argp.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include "data.h"
#include "io.h"
#include "ref.h"
#include "traj.h"

struct arguments {
    int n;              /* chain length */
    double temp;        /* temperature */
    char *fndata;       /* name of data file */
    char *fndump;       /* name of dump file */
    char *fntemp;       /* name of temp dump file */
};

const char *argp_program_version = "0.1";
const char *argp_program_bug_address = "Gustavo Rondina <rondina@gmail.com>";
static char doc[] = "qhe -- Quasi-harmonical approach to the entropy of macromolecules";

static struct argp_option options[] = {
    { "data",         'a', "FILE",  0, "LAMMPS data file",               0 },
    { "dump",         'u', "FILE",  0, "LAMMPS dump file",               0 },
    { "temp",         'e', "FILE",  0, "Temporary dump file",            0 },
    { "temperature",  't', "VALUE", 0, "Temperature",                    0 },
    { "chain-length", 'c', "VALUE", 0, "Length of Lennard-Jones chains", 0 },
    {  NULL,           0,   NULL,   0,  NULL,                            0 }
};

static error_t parse_opt(int key, char *arg, struct argp_state *state)
{
    struct arguments *arguments = state->input;

    switch (key) {
    case 'a':
        arguments->fndata = arg;
        break;
    case 'u':
        arguments->fndump = arg;
        break;
    case 'e':
        arguments->fntemp = arg;
        break;
    case 't':
        arguments->temp = atof(arg);
        break;
    case 'c':
        arguments->n = atoi(arg);
        break;
    case ARGP_KEY_ARG:
        argp_failure(state, 1, 0, "this program takes no arguments");
        break;
    default:
        return ARGP_ERR_UNKNOWN;
    }

    return 0;
}

static int check_arguments(struct arguments *arguments)
{
    int ret = 0;

    if (arguments->fndata == NULL) {
        fprintf(stderr, "Error: LAMMPS data file was not specified, see --help\n");
        ret = -1;
    }

    if (arguments->fndump == NULL) {
        fprintf(stderr, "Error: LAMMPS dump file was not specified, see --help\n");
        ret = -1;
    }

    if (arguments->fntemp == NULL) {
        fprintf(stderr, "Error: temporary dump file was not specified, see --help\n");
        ret = -1;
    }

    if (arguments->temp <= 0) {
        fprintf(stderr, "Error: temperature must be greater than zero\n");
        ret = -1;
    }

    if (arguments->n <= 0) {
        fprintf(stderr, "Error: chain length must be greater than zero\n");
        ret = -1;
    }

    return ret;
}

int main(int argc, char **argv)
{
    struct arguments arguments;
    arguments.n = -1;
    arguments.temp = -1;
    arguments.fndata = NULL;
    arguments.fndump = NULL;
    arguments.fntemp = NULL;

    struct argp argp = { options, parse_opt, NULL, doc, NULL, NULL, NULL };
    argp_parse(&argp, argc, argv, ARGP_IN_ORDER, 0, &arguments);

    if (check_arguments(&arguments) < 0)
        return EXIT_FAILURE;

    struct data data;
    init_data(&data);

    if (read_data(arguments.fndata, &data) < 0)
        return EXIT_FAILURE;

    /* Make sure info read from data file is valid */
    assert(data.natoms > 0);
    
    /* There should be at least n atoms (i.e., at least 1 molecule) */
    assert(arguments.n <= data.natoms);
    
    /* All molecules should have n atoms */
    assert((data.natoms % arguments.n) == 0);

    /* Fill other fields */
    data.temp = arguments.temp;
    data.molsize = arguments.n;
    data.nmols = data.natoms/data.molsize;

    /* Reference */
    struct molecule *refmols = init_reference(&data);

    /* 1st pass */
    parse_pass1(arguments.fndump, &data, refmols);

    /* 2nd pass */
    parse_pass2(arguments.fndump, arguments.fntemp, &data, refmols);

    /* Clean up */
    free_reference(refmols, &data);
    free_data(&data);

    return EXIT_SUCCESS;
}
