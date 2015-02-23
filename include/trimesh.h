#include <stdio.h>

#ifndef trimesh_h
#define trimesh_h

#include "Point2D.h"

int trimesh(int nn, int tdim, int nb, int nab, int *nbs, int ***bs, Point2D *pt, int tri[][3]);
int trimesh(int nn, int tdim, int nb, int nab, int *nbs, int ***bs, double *x, double *y, int tri[][3]);

#endif
