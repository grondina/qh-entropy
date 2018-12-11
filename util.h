#ifndef _UTIL_H
#define _UTIL_H

#include "data.h"

extern int cmpint(const void *p1, const void *p2);
extern double gyration(struct molecule *mol);
extern void kabsch(struct molecule *mol, struct molecule *ref);
extern void remove_com(struct molecule *mol);

#endif
