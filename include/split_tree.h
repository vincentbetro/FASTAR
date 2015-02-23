#include <stdio.h>
#include "geometry.h"
#include "List.h"

#define STCHUNK 10000

#ifndef split_tree_h
#define split_tree_h

class TREEVOXEL
{
 public:
  void print(FILE *outf)
  {
    fprintf(outf,"\nMother = %d",mom);
    fprintf(outf,"\nSplit  = %d",split);
    fprintf(outf,"\nCut = %d",cut);
    int j;
    
    fprintf(outf,"\nVoxel kids:");
    for (j=0; j < children.max; j++)
      fprintf(outf,"\nKid %d = %d",j,children.list[j]);
    
    fprintf(outf,"\nVoxel nodes:");
    for (j=0; j < 27; j++)
      if (nodes[j] >= 0)
        fprintf(outf,"\nNode %d = %d",j,nodes[j]);
        
    fprintf(outf,"\nVoxel cornerpoints:");
      for (j=0; j < 2; j++)
        fprintf(outf,"\nCornerpt %d = (%lf, %lf, %lf)",j,cornerpts[j][0],cornerpts[j][1],cornerpts[j][2]);
    fprintf(outf,"\n");
  }
  TREEVOXEL()
  {
    mom = -1;
    split = 0;
    cut = 0;
    for (int i = 0; i < 27; i++)
      nodes[i] = -1;
  }
  ~TREEVOXEL()
  {
    split = mom = cut = 0;
  }

  int split, mom, cut;
  Point cornerpts[2];
  List children;
  int nodes[27];
};

class split_tree
{
 public:
  split_tree()
  { nvoxels=voxel_dim=0;
    voxels = 0;
  }
  ~split_tree()
  {
    if (voxels != 0)
      free(voxels);
    nvoxels=voxel_dim=0;
  }
  
  void create_tree(geometry *geom, double &smn, double &smx, double ar, double mtol, double nspace);
  void cut_voxels(geometry *geom, double ctol, double smn, double smx, int *boxcut, int *vbd);
  void refine_voxel(int &new_voxel, int k, int code);
  int nabor_search(Point pt, int in, double smn, double smx);
  void set_nodes(int v, double smn, double smx, int &num_nodes, int *nums, Point *pts);
  void create_voxels(int in, double RMT[3][3], double smn, double smx, double ar, double hi[3], double lo[3], double mtol);
  int find_nabors_nodes(int nbr, double tol, Point pt);
  void set_nabors_nodes(int k, int nbr, double tol, Point pt);
  int eval_qualityCB(int v, int nbr, double del[3], int type[2], int &changes);
  void eval_quality41(int nbr, double del[3], int type[2]);
  void nabor_nodes(Point *pt, Point p0, Point p1);
  void voxel_search(geometry *geom, double hi[3], double lo[3], int in, double ctol, int f, int *vbd, int bd);
  void voxel_delete(geometry *geom, int in, double ctol, int f, int *vbd, int bd);
  int voxel_grad(List *nbrs0, List *nbrs1, int dir0, int dir1, double del[3], int i, int &changes3, double smn, double smx, Point pt1, Point pt2);
  void voxel_search_tree(Point tpt, int in, double ctol, List *out);

  int nvoxels,voxel_dim;
  TREEVOXEL *voxels;
};

#endif

