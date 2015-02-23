#include <stdio.h>

#ifndef Pmap_h
#define Pmap_h
class Pmap
{
 public:
  
  Pmap (int g = -1, int p = -1, int i = -1)
  {
    global = g;
    proc = p;
    index = i;
  }
  ~Pmap(){}
  void print(FILE *outf)
  {
    fprintf(outf,"\nGlobal= %d, Processor= %d, Index= %d",global,proc,index);
  }
  
  int global;
  int proc;
  int index;
};

inline bool operator==(Pmap a, Pmap b)
{
  return (a.proc == b.proc && a.index == b.index);
}

inline bool operator!=(Pmap a, Pmap b)
{
  return (a.proc != b.proc || a.index != b.index);
}

#endif
