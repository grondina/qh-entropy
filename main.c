#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <argp.h>
#include "io.h"

struct arguments {
    char *fndata;
    char *fndump;
};

const char *argp_program_version = "0.1";
const char *argp_program_bug_address = "Gustavo Rondina <rondina@gmail.com>";
static char doc[] = "qh-entropy -- Quasi-harmonical approach to the entropy of macromolecules";

static struct argp_option options[] = {
    { "data", 'a', "FILE", 0, "LAMMPS data file", 0 },
    { "dump", 'u', "FILE", 0, "LAMMPS dump file", 0 },
    { NULL,    0,  NULL,   0, NULL,               0 }
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

    return ret;
}

int main(int argc, char **argv)
{
    struct arguments arguments;

    arguments.fndata = NULL;
    arguments.fndump = NULL;

    struct argp argp = { options, parse_opt, NULL, doc, NULL, NULL, NULL };
    argp_parse(&argp, argc, argv, ARGP_IN_ORDER, 0, &arguments);

    if (check_arguments(&arguments) < 0)
        return EXIT_FAILURE;

    return EXIT_SUCCESS;
}
