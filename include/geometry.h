#include <stdio.h>
#include "Point.h"
#include "Vector.h"
#include "List.h"
#include "PList.h"
#include "Geom_Edge.h"
#include "Geom_Facet.h"
#include "Octree_Storage.h"

#ifndef mat_combo_h
#define mat_combo_h
class Mat_Combo
{
  public:
  Mat_Combo()
  { nm = 0; mat = 0; }
  ~Mat_Combo()
  { if (nm > 0) free(mat); nm = 0; mat = 0; }
  int Add_mat(int m)
  {
    if (mat == 0)
      mat = (int*)malloc(sizeof(int));
    else
      mat = (int*)realloc((void*)mat,(nm+1)*sizeof(int));
    mat[nm++] = m;

    return(nm);
  }

  int nm;
  int *mat;
};
#endif

#ifndef chain_h
#define chain_h
class Chain
{
  public:
  Chain()
  {  nn = nk = mat = 0; }
  ~Chain()
  { free(nodes); 
    if (nk > 0) free(knots);
    nn = nk = mat = 0;
  }

  void print(FILE *outf)
  {
    int n;
    fprintf(outf,"\nNumber of nodes = %d",nn);
    fprintf(outf,"\nNumber of knots = %d",nk);
    fprintf(outf,"\nMaterial = %d",mat);
    for (n=0;n < nn; n++)
      fprintf(outf,"\nNode %d = %d",n,nodes[n]);
    for (n=0;n < nk; n++)
      fprintf(outf,"\nKnot %d = %d",n,knots[n]);
  }

  int nn, nk;
  int *nodes, *knots;
  int mat;
};
#endif

#ifndef critical_point_h
#define critical_point_h
class Critical_Point
{
  public:
  Critical_Point()
  { nmat = 0; mat = 0; }
  ~Critical_Point()
  { free(mat); nmat = 0; }

  int nmat, index;
  int *mat;
};
#endif

#ifndef octree_list_h
#define octree_list_h
class Octree_List
{
 public:
  Octree_List(int l);
  ~Octree_List();

  void retrieve_list(Point p, Point tol, int lvl, List *list);
  int max_level(int l);
  double finest_size(Point p);

  Octree_List *kids[2][2][2];
  Octree_List *mom;
  Point lo, hi, mid;
  List* facets;
  int level;
};
#endif

#ifndef geometry_h
#define geometry_h
class geometry
{
 public:
  geometry()
  {
    //printf("\nConstructing geometry object.");
    ngb = 0;
    n_gnodes = 0;
    n_gfacets = 0;
    n_begin = 0;
    n_end = 0;
    layers = 0;
    g_space = 0;
    n_space = 0;
    g_vert = 0;
    blo = 0;
    bhi = 0;
    g_facet = 0;
    glist = 0;
    g_bname = 0;
    chain = 0;
    Octree_root = 0;
  }
  ~geometry()
  { 
    //printf("\nDestructing facet geometry object.\n");
    if (g_bname != 0 && ngb > 0)
    {
      for (int i=0; i < ngb; i++)
        free(g_bname[i]);
      free(g_bname);
    }
    if (g_vert != 0) free(g_vert);
    if (n_begin != 0) delete[] n_begin;
    if (n_end != 0) delete[] n_end;
    ngb=n_gnodes=n_gfacets=0;
    if (g_facet != 0) free(g_facet);
    if (g_space != 0) delete[] g_space;
    if (n_space != 0) delete[] n_space;
    if (layers != 0) delete[] layers;
    if (chain != 0) free(chain);
    if (blo != 0) free(blo);
    if (bhi != 0) free(bhi);
    if (glist != 0)
    {
      glist->Redimension(0);
      delete glist;
    }
    if (Octree_root != 0) delete Octree_root;
  }
  void read_geom(int ngf, char* fnames[]);
  int Identify_material_combos();
  int Intersect(Point pa, Point pb, Point &p);
  Point min_point();
  Point max_point();
  Point closest(Point p, int &b, Vector &nrm);
  Point closest_chain(Point p, int &b);
  void Plotfile(int mode);
  int In_Out(Point pt, double tol);
  
  int ngb, n_gnodes, n_gfacets, n_chain, Octree_max_level, n_combos, ncp;
  int *n_begin, *n_end, *layers;
  double *g_space, *n_space;
  char **g_bname;
  Point *g_vert, *blo, *bhi;
  Geom_Facet *g_facet;
  Chain *chain;
  Octree_Storage *Octree_root;
  PList *glist;
  Mat_Combo *mat_combo;
  Critical_Point *c_point;
};
#endif
