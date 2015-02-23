#include <stdio.h>
#include "List.h"

#ifndef Bar_n_h
#define Bar_n_h

class Bar_n
{
 public:

  void init(int n1, int n2)
  {
  nodes = 0;
  nodes = new List();
  nodes->Redimension(0);
  nodes->Add_To_List(n1);
  nodes->Add_To_List(n2);
  }
  void print(FILE *fp)
  {
    fprintf(fp,"\nEdge nodes:");
    int i;
    for (i=0; i < nodes->max; i++)
      fprintf(fp,"\nNode %d = %d",i,nodes->list[i]);
  }

  List *nodes;
};
#endif
