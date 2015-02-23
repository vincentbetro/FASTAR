#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include "geometry.h"
#include "split_tree.h"
#define MIN(x,y) ((x) <= (y) ? (x) : (y))
#define MAX(x,y) ((x) >= (y) ? (x) : (y))

void fastar_lib(geometry *geom, double smn, double smx, double ar, double ctol, double mtol, double nspace, int *boxcut, int *vbd, int *regions);		   
		   
