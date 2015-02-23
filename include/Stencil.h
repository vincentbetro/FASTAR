#include <stdio.h>
#include <stdlib.h>
#include "Point.h"
#include "List.h"

#ifndef Stencil_objects_h
#define Stencil_objects_h
class Stencil
{
 public:
  Stencil()
  {
    angle = 0;
    cnode = 0;
    nabor = 0;
    tri = 0;
    edge = 0;
    nt = 0;
    ne = 0;
  }
  ~Stencil()
  {
    if (angle != 0) delete[] angle;
    if (cnode != 0) delete[] cnode;
    if (nabor != 0) delete nabor;
    if (tri != 0) delete[] tri;
    if (edge != 0) delete[] edge;
    nt = 0;
    ne = 0;
  }
  void Initialize()
  {
    angle = 0;
    cnode = 0;
    nabor = 0;
    edge = 0;
    tri = 0;
    nt = 0;
    ne = 0;
  }
  void print(FILE *prnt)
  {
    int i, n;
    fprintf(prnt,"\nInput %d triangles",nt);
    for (i=0; i < nt; i++)
      fprintf(prnt,"\n  triangle %d contains nodes %d, %d, %d",i,tri[i][0],tri[i][1],tri[i][2]);
    fprintf(prnt,"\nInput %d edges",ne);
    for (i=0; i < ne; i++)
      fprintf(prnt,"\n  edge %d contains nodes %d, %d",i,edge[i][0],edge[i][1]);
    fprintf(prnt,"\nNabor list contains %d nodes",nabor->max);
    for (i=0; i < nabor->max; i++)
    {
      n = nabor->list[i];
      fprintf(prnt,"\n  neighbor %d = %d, cnode = (%g, %g, %g)",i,n,cnode[i][0],cnode[i][1],cnode[i][2]);
    }
  }

  int Create_Stencil(int ntri, int (*tri_conn)[3]);
  double Optimize_Stencil(int cflag, int niter, double threshold, FILE* prnt);
  double cost_function(int n, int cflag);
  void check_angles(int nn, double min_ang, double max_ang, const int mode);

  int nt, ne;
  int (*tri)[3];
  int (*edge)[2];
  double *angle;
  List *nabor;
  Point *cnode;
};
#endif
