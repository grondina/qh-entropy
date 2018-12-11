#ifndef _REF_H
#define _REF_H

#include "data.h"

extern struct molecule *init_reference(struct data *data);
extern void free_reference(struct molecule *refmols, struct data *data);
void print_reference(struct molecule *refmols, struct data *data);

#endif
