#ifndef _REF_H
#define _REF_H

#include "data.h"

extern struct molecule *refmols;

extern void init_reference(struct data *data);
extern void free_reference(struct data *data);

#endif
