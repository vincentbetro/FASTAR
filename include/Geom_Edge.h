#include <stdio.h>

#ifndef Bar_2_h
#define Bar_2_h
class Bar_2
{
 public:

  Bar_2(int n0 = -1, int n1 = -1) { nodes[0]=n0; nodes[1]=n1; }

  void print(FILE *outf)
  {
    fprintf(outf,"\nEdge nodes:");
    int i;
    for (i=0; i < 2; i++)
      fprintf(outf,"\nNode %d = %d",i,nodes[i]);
  }

  int nodes[2];
};
#endif

#ifndef Bar_3_h
#define Bar_3_h
class Bar_3
{
 public:

  Bar_3(int n0 = -1, int n1 = -1, int n2 = -1) { nodes[0]=n0; nodes[1]=n1; nodes[2]=n2;}

  void print(FILE *outf)
  {
    fprintf(outf,"\nEdge nodes:");
    int i;
    for (i=0; i < 3; i++)
      fprintf(outf,"\nNode %d = %d",i,nodes[i]);
  }

  int nodes[3];
};
#endif
