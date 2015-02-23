#include <stdio.h>
#include "Edge_n.h"

#ifndef face_h
#define face_h

class face
{
 public:

  void init()
  {
  edges = 0;
  bedges = 0;
  elements = 0;
  edges = new List();
  bedges = new List();
  elements = new List();
  edges->Redimension(0);
  bedges->Redimension(0);
  elements->Redimension(0);
  }
  void print(FILE *fp)
  {
    int i;
    fprintf(fp,"\nBoundary edges:");
    for (i=0; i < bedges->max; i++)
      fprintf(fp,"\nBoundary edge %d = %d",i,bedges->list[i]);
    fprintf(fp,"\nIncluded edges:");
    for (i=0; i < edges->max; i++)
      fprintf(fp,"\nEdge %d = %d",i,edges->list[i]);
    fprintf(fp,"\nIncluded in the following elements:");
    for (i=0; i < elements->max; i++)
      fprintf(fp,"\nElement %d = %d",i,elements->list[i]);
  }

  List *edges;
  List *bedges;
  List *elements;
};
#endif
