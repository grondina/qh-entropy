#include <stdio.h>
#include <stdlib.h>

int cmpint(const void *p1, const void *p2)
{
    int a = *((int *)p1);
    int b = *((int *)p2);
    return a - b;
}
