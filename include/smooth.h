#include <stdio.h>
#include "Pmap.h"
#include "Point.h"
#include "geometry.h"
#include "List.h"
#include "Bnormal.h"
#include "Stencil.h"
#include "Edge_n.h"
#include "face.h"

#ifndef Poly_mesh_objects_h
#define Poly_mesh_objects_h

#define EDGE_CHUNK 1000
#define FACE_CHUNK 1000
#define NODE_CHUNK 1000
#define ELEMENT_CHUNK 1000
#define TRI_CHUNK 1000

class POLY_ELEMENT
{
 public:
  POLY_ELEMENT() { nf = 0; f_n = 0; me = Pmap(-1,-1,-1); return;} // Constructor
  ~POLY_ELEMENT() // Destructor
  {
    if (nf > 0)
    {
      for (int f=0; f < nf; f++)
        delete f_n[f];
      delete[] f_n;
    }
    nf = 0;
    f_n = 0;
    me = Pmap(-1,-1,-1);
    return;
  }
  void Initialize() // for Resetting or intializing an array of elements.
  {
    nf = 0;
    f_n = 0;
    me = Pmap(-1,-1,-1);
    return;
  }
  void add_face(List *node_list) // add polygon as a face
  {
    if (node_list->max == 0)
      return;
    f_n = (List**)realloc((void*)f_n,(nf+1)*sizeof(List*));
    f_n[nf] = new List();
    *f_n[nf] = *node_list;
    nf++;
    return;
  }

  int nf;     // number of faces in polyhedral element (1 is a polygon)
  List **f_n; // list of nodes per face
  Pmap me;    // processor map
};

class NODE
{
 public:
  NODE()
  {
    stencil.Initialize();
    me = Pmap(-1,-1,-1);
    return;
  }
  ~NODE()
  {
    stencil.Initialize();
    me = Pmap(-1,-1,-1);
    return;
  }
  void print(FILE *prnt)
  {
    vert.print(prnt);
    me.print(prnt);
    stencil.print(prnt);
  }

  Point vert;      // Coordinate triplet
  Stencil stencil; // Virtual Control Volume
  Pmap me;         // Parallel map info
};

class BOUNDARY_MESH
{
 public:
  BOUNDARY_MESH() { elist=0; name=0; return; }
  ~BOUNDARY_MESH()
  {
    if (elist > 0) delete elist;
    if (name > 0) delete[] name;
    return;
  }

  List *elist;
  char *name;
};

class Vmesh_obj
{
 public:
  Vmesh_obj()
  { nfaces=nedges=edge_dim=face_dim=0;
    faces = 0;
    edges = 0;
  }
  ~Vmesh_obj()
  {
    if (faces != 0)
      free(faces);
    if (edges != 0)
      free(edges);
    nfaces=nedges=edge_dim=face_dim=0;
  }

  int nedges, nfaces, face_dim, edge_dim;
  face *faces;
  Bar_n *edges;
};

class POLYMESH
{
 public:
  POLYMESH()
  {
    nn=nelems=nb=nghost=element_dim=node_dim=0;
    node=0;
    element=0;
    boundary=0;
    pmap=0;
    vmesh = 0;
    return;
  }
  ~POLYMESH()
  {
    //if (node != 0) free(node);
    //if (element != 0) free(element);
    //if (boundary != 0) free(boundary);
    if (vmesh != 0) free(vmesh);
    if (pmap != 0)
      {
	  for (int i = 0; i < nn; i++)
	    free(pmap[i]);
	  free(pmap);
      }
    nn=nelems=nb=element_dim=nghost=node_dim=0;
    node=0;
    element=0;
    boundary=0;
    pmap=0;
    vmesh = 0;
    return;
  }
  Point min_point();
  Point max_point();
  int adjust_VCV_angles(int mode, double sa_min, double sa_max, double factor);
  void boundary_condition(int n, double u[], double v[], double w[]);
  int check_volumes();
  int check_Jacobians();
  bool control_volume_io(char sname[],int mode);
  int create_compress_row_storage(int **ia, int **iau, int **ja);
  void create_control_volumes();
  void create_element_maps();
  void create_viscous_layers(int v_layers, geometry *geom);
  int create_floating_nodes(int wsmoo, int ftag[], Bnormal **bf_nodes);
  void create_node_to_node(List **nhash);
  int critical_points(geometry *geom, int **cpn);
  void evaluate_control_volumes();
  void gg_gradients(int n, double f[], double &fx, double &fy, double &fz);
  void mesh_stats();
  int read_mesh(char sname[]);
  int write_mesh(char sname[]);
  int smooth_io(int mode, geometry *geom, char sname[]);
  void exchange_node_double(double *array);
  void exchange_node_int(int *array);
  void Winslow(geometry *geom, int osmoo, int nsmoo, int nsearch, int cnvrg, int mvbnd,
               int v_layers, double geo1, double geo2, double geo_min, double geo_max,
               int VCV_angle, double VCV_amin, double VCV_amax, double VCV_edge);
  void build_Winslow_matrix(int tag[], double mat[][3][3], double u[], double v[], double w[],
                                    int ia[], int iau[], int ja[]);
  void stencil_plotfile(int nd, Point cg, const char *prefix);
  void set_up_tri(int &tnt, List *tnode, int *map);
  void set_up_tri_elem(int &tnt, List *tnode, int *map);

  int nn, nelems, nb, nghost;
  int **pmap;
  NODE *node;
  POLY_ELEMENT *element;
  BOUNDARY_MESH *boundary;
  Vmesh_obj *vmesh;
  int element_dim, node_dim; //dimension of elem/node arrays
};
#endif
