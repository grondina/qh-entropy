#include <argp.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>

#include "data.h"
#include "io.h"
#include "traj.h"
#include "util.h"

const char *argp_program_version = "0.1";
const char *argp_program_bug_address = "Gustavo Rondina <rondina@gmail.com>";
static char doc[] = "qhe -- Quasi-harmonical approach to the entropy of macromolecules";

static struct argp_option options[] = {
    { "data",         'a', "FILE",  0, "LAMMPS data file",               0 },
    { "dump",         'u', "FILE",  0, "LAMMPS dump file",               0 },
    { "temp",         'e', "FILE",  0, "Temporary dump file",            0 },
    { "temperature",  't', "VALUE", 0, "Temperature",                    0 },
    { "chain-length", 'c', "VALUE", 0, "Length of Lennard-Jones chains", 0 },
    { "boltzmann",    'k', "VALUE", 0, "Boltzmann constant",             0 },
    { "planck",       'h', "VALUE", 0, "Planck's constant",              0 },
    { "chain-length", 'c', "VALUE", 0, "Length of Lennard-Jones chains", 0 },
    { "out",          'o', "FILE",  0, "File to write final entropies",  0 },
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
    case 'o':
        arguments->fnsave = arg;
        break;
    case 't':
        arguments->temp = atof(arg);
        break;
    case 'c':
        arguments->n = atoi(arg);
        break;
    case 'k':
        arguments->kB = atof(arg);
        break;
    case 'h':
        arguments->h = atof(arg);
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

    if (arguments->fnsave == NULL) {
        fprintf(stderr, "Error: output file was not specified, see --help\n");
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

    if (arguments->kB < 0) {
        fprintf(stderr, "Error: Boltzmann's constant must be greater than zero\n");
        ret = -1;
    }

    if (arguments->h < 0) {
        fprintf(stderr, "Error: Planck's constant must be greater than zero\n");
        ret = -1;
    }

    return ret;
}

static void init_arguments(struct arguments *arguments)
{
    arguments->n = -1;
    arguments->h = -1;
    arguments->kB = -1;
    arguments->temp = -1;
    arguments->fndata = NULL;
    arguments->fndump = NULL;
    arguments->fntemp = NULL;
    arguments->fnsave = NULL;
}

int main(int argc, char **argv)
{
    struct arguments arguments;
    init_arguments(&arguments);

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
    struct molecule *refmols = init_molecule_array(&data);

    /* Average */
    struct molecule *avemols = init_molecule_array(&data);

    /* 1st pass */
    parse_pass1(arguments.fndump, &data, refmols);

    /* 2nd pass */
    parse_pass2(arguments.fndump, arguments.fntemp, &data, refmols, avemols);

    /* 3nd pass */
    int n = 3 * data.molsize;
    double (*sigma)[n][n] = malloc(sizeof(double[data.nmols][n][n]));;
    double (    *M)[n][n] = malloc(sizeof(double[data.nmols][n][n]));
    assert(sigma != NULL);
    assert(    M != NULL);
    parse_pass3(arguments.fntemp, &data, avemols, n, sigma, M);

    /* Calculate entropy */
    double *S = entropy(n, &data, &arguments, sigma, M);
    
    /* Print & save entropy */
    write_entropy(arguments.fnsave, &data, S);

    /* Clean up */
    free(sigma);
    free(M);
    free(S);
    free_molecule_array(refmols, &data);
    free_molecule_array(avemols, &data);
    free_data(&data);

    return EXIT_SUCCESS;
}
