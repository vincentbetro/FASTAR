//#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "geometry.h"
#include "List.h"
#include "Linked_List.h"
#include "Vector.h"
#include "Util.h"
#include "PList.h"
#include "Octree_Storage.h"

extern int my_rank;
extern int num_procs;
extern FILE *out_f;

#define FACET_TYPE 1
#define TRI_TYPE   2

void geometry::read_geom(int ngf, char* gnames[])
{
  fprintf(out_f,"\nNumber of geometry files = %d",ngf);

  int mat_max = 0;
  const int bdim = 132;
  char buff[bdim], fname[bdim];
  FILE* g_io;
  int c, f, g, i, j, k, m, mat, n, nf, nfc, npts, n0, n1, n2, n3, s;
  int ftype;
  int  **itmp;

  ngb = n_gnodes = 0;

  g_facet = 0;

  // Read each geometry file
  for (i=0; i < ngf; i++)
  {
    fprintf(out_f,"\nGeometry file %d = %s",i+1,gnames[i]);
    if ((g_io=fopen(gnames[i],"r")) == 0)
    {
      fprintf(stderr,"\nCouldn't open file %s",gnames[i]);
      fflush(stderr);
      //MPI_Abort(MPI_COMM_WORLD,0);
      exit(0);
    }

    double x, y, z;
    // determine geometry file type by suffix
    ftype = 0;
    if (strstr(gnames[i],".facet") != NULL)
      ftype=FACET_TYPE;
    if (strstr(gnames[i],".tri") != NULL)
      ftype=TRI_TYPE;
    switch(ftype)
    {
      case (FACET_TYPE):
        fgets(buff,bdim,g_io);
        fgets(buff,bdim,g_io);
        fgets(buff,bdim,g_io);
        fgets(buff,bdim,g_io);
        fgets(buff,bdim,g_io);
        sscanf(buff,"%d",&npts);
        fprintf(out_f,"\n Number of points = %d",npts);
        if (g_vert != 0)
          g_vert=(Point*)realloc((void*)g_vert,
                                  (n_gnodes+npts)*sizeof(Point));
        else
          g_vert=(Point*)malloc((n_gnodes+npts)*sizeof(Point));
        for (j=0; j < npts; j++)
        {
          fgets(buff,bdim,g_io);
          sscanf(buff,"%lf %lf %lf", &x, &y, &z);
          g_vert[n_gnodes+j] = Point(x,y,z);
        }
        fgets(buff,bdim,g_io);
        sscanf(buff,"%d",&nf);
        fprintf(out_f,"\n Number of faces = %d",nf);
        for (f=0; f < nf; f++)
        {
          fgets(buff,bdim,g_io);
          fgets(buff,bdim,g_io);
          sscanf(buff,"%d %d",&j,&k);
          nfc = (j == 3 && k != 3) ? k : j;
          fprintf(out_f,"\n Number of facets = %d",nfc);
          if (g_facet != 0)
            g_facet=(Geom_Facet*)realloc((void*)g_facet,
                                      (n_gfacets+nfc)*sizeof(Geom_Facet));
          else
            g_facet=(Geom_Facet*)malloc((n_gfacets+nfc)*sizeof(Geom_Facet));
          for (j=0; j < nfc; j++)
          {
            fgets(buff,bdim,g_io);
            sscanf(buff,"%d %d %d %d",&n1,&n2,&n3,&mat);
            mat_max = MAX(mat_max,mat);
            g_facet[n_gfacets+j].nodes[0] = n1-1+n_gnodes;
            g_facet[n_gfacets+j].nodes[1] = n2-1+n_gnodes;
            g_facet[n_gfacets+j].nodes[2] = n3-1+n_gnodes;
            g_facet[n_gfacets+j].mat = mat;
          }
          n_gfacets += nfc;
        }
        n_gnodes += npts;
        break;
      case (TRI_TYPE):
        fgets(buff,bdim,g_io);
        sscanf(buff,"%d %d",&npts,&nf);
        fprintf(out_f,"\n Number of points = %d",npts);
        fprintf(out_f,"\n Number of triangles = %d",nf);
        if (g_vert != 0)
          g_vert=(Point*)realloc((void*)g_vert,
                                  (n_gnodes+npts)*sizeof(Point));
        else
          g_vert=(Point*)malloc((n_gnodes+npts)*sizeof(Point));
        if (g_facet != 0)
          g_facet=(Geom_Facet*)realloc((void*)g_facet,
                                    (n_gfacets+nf)*sizeof(Geom_Facet));
        else
          g_facet=(Geom_Facet*)malloc((n_gfacets+nf)*sizeof(Geom_Facet));
        for (j=0; j < npts; j++)
        {
          fgets(buff,bdim,g_io);
          sscanf(buff,"%lf %lf %lf", &x, &y, &z);
          g_vert[n_gnodes+j] = Point(x,y,z);
        }
        for (j=0; j < nf; j++)
        {
          fgets(buff,bdim,g_io);
          sscanf(buff,"%d %d %d",&n1,&n2,&n3);
          g_facet[n_gfacets+j].nodes[0] = n1-1+n_gnodes;
          g_facet[n_gfacets+j].nodes[1] = n2-1+n_gnodes;
          g_facet[n_gfacets+j].nodes[2] = n3-1+n_gnodes;
        }
        for (j=0; j < nf; j++)
        {
          fgets(buff,bdim,g_io);
          sscanf(buff,"%d",&mat);
          mat_max = MAX(mat_max,mat);
          g_facet[n_gfacets+j].mat = mat;
        }
        n_gfacets += nf;
        n_gnodes += npts;
        break;
      default:
        break;
    }

    fclose(g_io);
  }
  fflush(out_f);

  Vector v1, v2, norm;
  int *map;
  // Eliminate duplicate nodes
  if (ngf > 1)
  {
    map = new int[n_gnodes];

    // Identify unused nodes
    fprintf(out_f,"\nChecking for unused nodes");
    for (n=0; n < n_gnodes; n++)
      map[n]= -1;
    for (f=0; f < n_gfacets; f++)
      for (i=0; i < 3; i++)
      {
        n = g_facet[f].nodes[i];
        map[n] = 1;
      }
    for (j=n=0; n < n_gnodes; n++)
      if (map[n]== -1)
        j++;
    if (j > 0)
      fprintf(out_f,"\nNumber of unused nodes = %d",j);

    // Identify duplicate nodes
    fprintf(out_f,"\nChecking for duplicate nodes");
    fprintf(out_f,"\n Each '.' equals 100 nodes\n");
    for (n=0; n < n_gnodes; n++)
      if (map[n] > 0)
        map[n]=n;
    for (j=n_gnodes-1; j > 1; j--)
    {
      if (j % 100 == 0)
        fprintf(out_f,".");
      if (j % 1000 == 0)
        fprintf(out_f," %d\n",j);
      fflush(stdout);
      for (i=0; i < j && map[j]==j; i++)
      {
        if (map[i] < 0)
          continue;
        if (distance(g_vert[i],g_vert[j]) < 1.0e-20)
          map[j] = i;
      }
    }
    // Update facet connectivity
    for (f=0; f < n_gfacets; f++)
      for (i=0; i < 3; i++)
      {
        n = g_facet[f].nodes[i];
        g_facet[f].nodes[i] = map[n];
      }
    // compress node list
    for (j=i=0; i < n_gnodes; i++)
    {
      if (map[i] == i)
      {
        if (j != i)
          g_vert[j] = g_vert[i];
        map[i] = j++;
      } else
        map[i] = -1;
    }
    if (j < n_gnodes)
    {
      fprintf(out_f,"\nNumber of geometry nodes deleted = %d",n_gnodes-j);
      g_vert = (Point*)realloc((void*)g_vert,j*sizeof(Point));
      n_gnodes = j;
      // Update facet connectivity
      for (f=0; f < n_gfacets; f++)
        for (i=0; i < 3; i++)
        {
          n = g_facet[f].nodes[i];
          g_facet[f].nodes[i] = map[n];
        }
    }
    delete[] map;

    map = new int[n_gfacets];
    // Eliminate zero area facets
    for (f=0; f < n_gfacets; f++)
    {
      n0 = g_facet[f].nodes[0];
      n1 = g_facet[f].nodes[1];
      n2 = g_facet[f].nodes[2];
      v1 = Vector(g_vert[n0],g_vert[n1]);
      v2 = Vector(g_vert[n0],g_vert[n2]);
      norm = v1 % v2;
      if (norm.magnitude() > 1.0e-20)
        map[f] = 1;
      else
        map[f] = 0;
    }
    for (j=i=0; i < n_gfacets; i++)
    {
      if (map[i] == 1)
      {
        if (j != i)
          g_facet[j] = g_facet[i];
        j++;
      }
    }
    delete[] map;

    if (j < n_gfacets)
    {
      fprintf(out_f,"\nNumber of geometry facets deleted = %d",n_gfacets-j);
      g_facet=(Geom_Facet*)realloc((void*)g_facet,j*sizeof(Geom_Facet));
      n_gfacets=j;
    }
  }
  fflush(out_f);

  // Count facets per material
  int* mat_map = new int[mat_max+1];
  for (i=0; i <= mat_max; i++)
    mat_map[i] = 0;
  for (i=0; i < n_gfacets; i++)
    mat_map[g_facet[i].mat]++;
  fprintf(out_f,"\nMaximum material ID = %d",mat_max);
  fflush(out_f);

  // Bodies are stored starting from 1
  // Virtual bodies will use location 0
  ngb = 0;
  for (i=0; i <= mat_max; i++)
    if (mat_map[i] > 0)
    {
      fprintf(out_f,"\nMaterial %d has %d facets.",i,mat_map[i]);
      if (i > 0)
        ngb = ngb + 1;
    }
  fflush(out_f);

  // create temporary to hold data while facets are reordered based on material
  itmp = new int*[5];
  for (i=0; i < 5; i++)
  {
    itmp[i] = new int[n_gfacets];
    if (!itmp[i])
    {
      fprintf(out_f,"\nREAD_GEOM: Unable to allocate space for itmp[%d]",i);
      fflush(out_f);
      //MPI_Abort(MPI_COMM_WORLD,0);
      exit(0);
    }
  }

  // Create body limits points
  blo = (Point*)malloc((mat_max+1)*sizeof(Point));
  bhi = (Point*)malloc((mat_max+1)*sizeof(Point));
  for (i=0; i <= mat_max; i++)
  {
    blo[i][0] = blo[i][1] = blo[i][2] =  1.0e20;
    bhi[i][0] = bhi[i][1] = bhi[i][2] = -1.0e20;
  }

  // Store facets in material order
  n_begin = new int[ngb+1];
  n_end = new int[ngb+1];
  g_bname = (char**)malloc((ngb+1)*sizeof(char*));
  for (i=0; i <= ngb; i++)
    g_bname[i] = (char*)malloc(33*sizeof(char));
  ngb = 0;
  nfc = 0;
  for (i=0; i <= mat_max; i++)
  {
    if (i==0 || mat_map[i] > 0)
    {
      if (i > 0) ngb++;
      n_begin[ngb] = nfc;
      k = 0;
      for (j=0; j < n_gfacets; j++)
        if (g_facet[j].mat == i)
        {
          k++;
          itmp[1][nfc] = g_facet[j].nodes[0];
          itmp[2][nfc] = g_facet[j].nodes[1];
          itmp[3][nfc] = g_facet[j].nodes[2];
          itmp[4][nfc] = ngb;
          itmp[0][nfc] = nfc;
          nfc++;
          int g0 = g_facet[j].nodes[0];
          int g1 = g_facet[j].nodes[1];
          int g2 = g_facet[j].nodes[2];
          Point p0 = g_vert[g0];
          Point p1 = g_vert[g1];
          Point p2 = g_vert[g2];
          blo[ngb][0] = MIN(blo[ngb][0],MIN(MIN(p0[0],p1[0]),p2[0]));
          blo[ngb][1] = MIN(blo[ngb][1],MIN(MIN(p0[1],p1[1]),p2[1]));
          blo[ngb][2] = MIN(blo[ngb][2],MIN(MIN(p0[2],p1[2]),p2[2]));
          bhi[ngb][0] = MAX(bhi[ngb][0],MAX(MAX(p0[0],p1[0]),p2[0]));
          bhi[ngb][1] = MAX(bhi[ngb][1],MAX(MAX(p0[1],p1[1]),p2[1]));
          bhi[ngb][2] = MAX(bhi[ngb][2],MAX(MAX(p0[2],p1[2]),p2[2]));
        }
      n_end[ngb] = n_begin[ngb] + k;
      fprintf(out_f,"\nBody %d has %d facets.",ngb,k);
    }
  }
  fprintf(out_f,"\nTotal number of facets = %d.",n_gfacets);
  fflush(out_f);

  delete[] mat_map;

  blo = (Point*)realloc((void*)blo,(ngb+1)*sizeof(Point));
  bhi = (Point*)realloc((void*)bhi,(ngb+1)*sizeof(Point));

  // Restore facets to structure and initialize facet spacing
  for (i=0; i < n_gfacets; i++)
  {
    j = itmp[0][i];
    n1 = itmp[1][i];
    n2 = itmp[2][i];
    n3 = itmp[3][i];
    g_facet[j].nodes[0] = n1;
    g_facet[j].nodes[1] = n2;
    g_facet[j].nodes[2] = n3;
    g_facet[j].mat = itmp[4][i];
  }

  // Free up temporary array
  for (i=0; i < 5; i++)
    delete[] itmp[i];
  delete[] itmp;

  // write new reduced facet file
  if (my_rank == 0 && ngf > 1)
  {
    if ((g_io=fopen("P_HUGG_reduced.facet","w")) == 0)
    {
      fprintf(stderr,"\nCouldn't open file P_HUGG_reduced.facet");
      fflush(stderr);
      //MPI_Abort(MPI_COMM_WORLD,0);
      exit(0);
    }
    fprintf(g_io,"Reduced Facet file written from P_HUGG\n");
    fprintf(g_io,"    1\n");
    fprintf(g_io,"Volume\n");
    fprintf(g_io,"0\n");
    fprintf(g_io,"%d\n",n_gnodes);
    for (n=0; n < n_gnodes; n++)
      fprintf(g_io,"%14.7e %14.7e %14.7e\n",g_vert[n][0],g_vert[n][1],g_vert[n][2]);
    if (n_end[0]-n_begin[0] == 0)
      fprintf(g_io,"%d\n",ngb);
    else
      fprintf(g_io,"%d\n",ngb+1);
    for (k=0; k <= ngb; k++)
    {
      if (k == 0 && n_end[0]-n_begin[0] == 0)
        continue;
      i = n_begin[k];
      j = g_facet[i].mat;
      i = n_end[k] - n_begin[k];
      if (i == 0) continue;
      fprintf(g_io,"Material %d\n",j);
      fprintf(g_io,"      3%7d      0      0      0      0      0\n",i);
      for (j=n_begin[k]; j < n_end[k]; j++)
      {
        n1 = g_facet[j].nodes[0]+1;
        n2 = g_facet[j].nodes[1]+1;
        n3 = g_facet[j].nodes[2]+1;
        m = g_facet[j].mat;
        fprintf(g_io,"%7d %7d %7d %7d\n",n1,n2,n3,m);
      }
    }
    fclose(g_io);
  }

  // Store chains forming the boundary of each boundary part

  // hash table stores faces per node
  List **fhash;
  fhash = (List**)malloc(n_gnodes*sizeof(List*));
  if (!fhash)
  {
    fprintf(out_f,"\nREAD_GEOM: unable to allocate memory for fhash!");
    fflush(out_f);
    //MPI_Abort(MPI_COMM_WORLD,0);
    exit(0);
  }
  for (n=0; n < n_gnodes; n++)
  {
    fhash[n] = new List();
    if (!fhash[n])
    {
      fprintf(out_f,"\nREAD_GEOM: Unable to allocate space for fhash[%d]",n);
      fflush(out_f);
      //MPI_Abort(MPI_COMM_WORLD,0);
      exit(0);
    }
  }
  for (g=0; g < n_gfacets; g++)
  {
    n0 = g_facet[g].nodes[0];
    n1 = g_facet[g].nodes[1];
    n2 = g_facet[g].nodes[2];
    fhash[n0]->Check_List(g);
    fhash[n1]->Check_List(g);
    fhash[n2]->Check_List(g);
  }

  // determine neighboring faces using hash table
  itmp = (int**)malloc(3*sizeof(int*));
  for (i=0; i < 3; i++)
  {
    itmp[i] = (int*)malloc(n_gfacets*sizeof(int));
    if (!itmp[i])
    {
      fprintf(out_f,"\nREAD_GEOM: Unable to allocate space for itmp[%d]",i);
      fflush(out_f);
      //MPI_Abort(MPI_COMM_WORLD,0);
      exit(0);
    }
  }
  for (g=0; g < n_gfacets; g++)
  {
    for (s=0; s < 3; s++)
    {
      switch (s)
      {
        case 0: 
          n0 = g_facet[g].nodes[0];
          n1 = g_facet[g].nodes[1];
          break;
        case 1: 
          n0 = g_facet[g].nodes[1];
          n1 = g_facet[g].nodes[2];
          break;
        case 2: 
          n0 = g_facet[g].nodes[2];
          n1 = g_facet[g].nodes[0];
          break;
      }
      k=-1;
      for (i=0; i < fhash[n0]->max && k == -1; i++)
      {
        m=fhash[n0]->list[i];
        if (m == g)
          continue;
        for (j=0; j < fhash[n1]->max && k == -1; j++)
          if (m == fhash[n1]->list[j])
            k = m;
      }
      itmp[s][g] = k;
    }
  }
  
  for (n=0; n < n_gnodes; n++)
    delete fhash[n];
  free(fhash);

  n_chain = 0;
  chain = 0;
  for (g=1; g <= ngb; g++)
  {
    int n_e = 0;
    Bar_2 *s_edge;
    s_edge = 0;
    for (i=n_begin[g]; i < n_end[g]; i++)
    {
      for (s=0; s < 3; s++)
      {
        switch (s)
        {
          case 0: 
            n0 = g_facet[i].nodes[0];
            n1 = g_facet[i].nodes[1];
            break;
          case 1: 
            n0 = g_facet[i].nodes[1];
            n1 = g_facet[i].nodes[2];
            break;
          case 2: 
            n0 = g_facet[i].nodes[2];
            n1 = g_facet[i].nodes[0];
            break;
        }
        m = itmp[s][i];
        if (g_facet[i].mat != g_facet[m].mat)
        {
          if (s_edge != 0)
            s_edge = (Bar_2*)realloc((void*)s_edge,(n_e+1)*sizeof(Bar_2));
          else
            s_edge = (Bar_2*)malloc((n_e+1)*sizeof(Bar_2));
          s_edge[n_e].nodes[0] = n0;
          s_edge[n_e].nodes[1] = n1;
          n_e++;
        }
      }
    }
    
    //
    // now store edges in contiguous strings
    //

    int flag = 0;

    while(n_e > 0)
    {
      if (flag == 0) // start a new string
      {
        if (chain != 0)
          chain = (Chain*)realloc((void*)chain,(n_chain+1)*sizeof(Chain));
        else
          chain = (Chain*)malloc((n_chain+1)*sizeof(Chain));
        chain[n_chain].nodes = (int*)malloc(2*sizeof(int));
        chain[n_chain].mat = g;
        chain[n_chain].nodes[0] = s_edge[n_e-1].nodes[0];
        chain[n_chain].nodes[1] = s_edge[n_e-1].nodes[1];
        chain[n_chain].nn = 2;
        n_e--;
        flag=1;
        if (n_e == 0)
          n_chain++;
      } else // attempt to add to either end of current string
      {
        int nn = chain[n_chain].nn-1;
        j = n0 = n1 = -1;
        for (i=0; i < n_e && j == -1; i++)
        {
          if (s_edge[i].nodes[0] == chain[n_chain].nodes[0])
          {
            n0 = s_edge[i].nodes[1];
            j = i;
          }
          if (s_edge[i].nodes[1] == chain[n_chain].nodes[0])
          {
            n0 = s_edge[i].nodes[0];
            j = i;
          }
          if (s_edge[i].nodes[0] == chain[n_chain].nodes[nn])
          {
            n1 = s_edge[i].nodes[1];
            j = i;
          }
          if (s_edge[i].nodes[1] == chain[n_chain].nodes[nn])
          {
            n1 = s_edge[i].nodes[0];
            j = i;
          }
        }
        if (n0 != -1) // Add to beginning of string
        {
          s_edge[j] = s_edge[--n_e];
          chain[n_chain].nodes = (int*)realloc((void*)chain[n_chain].nodes,
                                                 (nn+2)*sizeof(int));
          for (j=nn+1; j > 0; j--) // shift current node list up
            chain[n_chain].nodes[j] = chain[n_chain].nodes[j-1];
          chain[n_chain].nodes[0] = n0;
          chain[n_chain].nn++;
          if (n_e == 0)
            n_chain++;
        } else if (n1 != -1) // Add to end of string
        {
          s_edge[j] = s_edge[--n_e];
          chain[n_chain].nodes = (int*)realloc((void*)chain[n_chain].nodes,
                                                 (nn+2)*sizeof(int));
          chain[n_chain].nodes[nn+1] = n1;
          chain[n_chain].nn++;
          if (n_e == 0)
            n_chain++;
        } else  // Terminate current string
        {
          n_chain++;
          flag = 0;
        }
      }
    }
    free(s_edge);
  }
  fprintf(out_f,"\nNumber of chains = %d",n_chain);
  fflush(out_f);

  // Free up temporary array
  for (i=0; i < 3; i++)
    free(itmp[i]);
  free(itmp);

  // Determine knot locations
  for (c=0; c < n_chain; c++)
  {
    fprintf(out_f,"\nChain %d has body # = %d",c+1,chain[c].mat);
    chain[c].knots = 0;
    chain[c].nk = 0;
    int nn = chain[c].nn;
    if (chain[c].nodes[0] != chain[c].nodes[nn-1])
    {
      fprintf(out_f,"\ngeometry: Beginning and ending points of chain %d do not match!",c);
      fflush(out_f);
      //MPI_Abort(MPI_COMM_WORLD,0);
      exit(0);
    }
    for (n=0; n < nn-1; n++)
    {
      if (n==0)
        i = chain[c].nodes[nn-2];
      else
        i = chain[c].nodes[n-1];
      j = chain[c].nodes[n];
      k = chain[c].nodes[n+1];
      v1 = Vector(g_vert[i],g_vert[j]);
      v1.normalize();
      v2 = Vector(g_vert[j],g_vert[k]);
      v2.normalize();
      if (v1*v2 < 0.95)
      {
        int nk = chain[c].nk;
        if (chain[c].knots != 0)
          chain[c].knots = (int*)realloc((void*)chain[c].knots,(nk+1)*sizeof(int));
        else
          chain[c].knots = (int*)malloc((nk+1)*sizeof(int));
        chain[c].knots[nk] = j;
        chain[c].nk++;
      }
    }
    if (chain[c].nk > 0)
      fprintf(out_f,"\nChain %d has %d knots.",c+1,chain[c].nk);
  }
  fflush(out_f);

  // identify critical points in geometry (points with more than 2 materials)
  itmp = (int**)malloc(2*sizeof(int*));
  for (i=0; i < 2; i++)
  {
    itmp[i] = (int*)malloc(n_gnodes*sizeof(int));
    if (!itmp[i])
    {
      fprintf(out_f,"\nREAD_GEOM: Unable to allocate space for itmp[%d]",i);
      fflush(out_f);
      //MPI_Abort(MPI_COMM_WORLD,0);
      exit(0);
    }
  }

  for (n=0; n < n_gnodes; n++)
    itmp[0][n] = 0;
  for (i=1; i <= ngb; i++)
  {
    for (n=0; n < n_gnodes; n++)
      itmp[1][n] = 0;
    for (j=n_begin[i]; j < n_end[i]; j++)
    {
      n1 = g_facet[j].nodes[0];
      n2 = g_facet[j].nodes[1];
      n3 = g_facet[j].nodes[2];
      itmp[1][n1] = 1;
      itmp[1][n2] = 1;
      itmp[1][n3] = 1;
    }
    for (n=0; n < n_gnodes; n++)
      itmp[0][n] += itmp[1][n];
  }
  ncp = 0;
  for (n=0; n < n_gnodes; n++)
    if (itmp[0][n] > 2)
      ncp++;

  if (ncp > 0)
  {
    c_point = (Critical_Point*)malloc(ncp*sizeof(Critical_Point));
    ncp = 0;
    for (n=0; n < n_gnodes; n++)
    {
      if (itmp[0][n] > 2)
      {
        c_point[ncp].mat = (int*)malloc(itmp[0][n]*sizeof(int));
        c_point[ncp].nmat = 0;
        c_point[ncp].index = n;
        itmp[1][n] = ncp;
        ncp++;
      } else
        itmp[1][n] = -1;
    }
    for (i=1; i <= ngb; i++)
    {
      for (j=n_begin[i]; j < n_end[i]; j++)
      {
        for (k=0; k < 3; k++)
        {
          n = g_facet[j].nodes[k];
          n1 = itmp[1][n];
          if (n1 >= 0)
          {
            n2=-1;
            for (n3=0; n3 < c_point[n1].nmat && n2 < 0; n3++)
            {
              if (c_point[n1].mat[n3] == i)
                n2 = n3;
            }
            if (n2 < 0)
            {
              c_point[n1].mat[c_point[n1].nmat] = i;
              c_point[n1].nmat++;
            }
          }
        }
      }
    }
  }

  for (i=0; i < ncp; i++)
  {
    Point p = g_vert[c_point[i].index];
    fprintf(out_f,"\n\nCritical point %d is located at (%g, %g, %g)",i+1,p[0],p[1],p[2]);
    for (j=0; j < c_point[i].nmat; j++)
      fprintf(out_f,"\n  Boundary %d is %d",j+1,c_point[i].mat[j]);
  }
  fflush(out_f);

  // Free up temporary array
  for (i=0; i < 2; i++)
    free(itmp[i]);
  free(itmp);

  // compute area sums for each body and total
  double tx, ty, tz, bx, by, bz, twet;
  tx = ty = tz = twet = 0.0;
  double *wet, *perimeter;
  wet = new double[ngb+1];
  perimeter = new double[ngb+1];

  fprintf(out_f,"\n\nArea Summations");
  fprintf(out_f,"\nBoundary #     Area(x)        Area(y)        Area(z)        Wetted");
  fprintf(out_f,"\n---------- -------------- -------------- -------------- --------------");
  for (i=0; i <= ngb; i++)
  {
    bx = by = bz = wet[i] = perimeter[i] = 0.0;
    for (j=n_begin[i]; j < n_end[i]; j++)
    {
      n1 = g_facet[j].nodes[0];
      n2 = g_facet[j].nodes[1];
      n3 = g_facet[j].nodes[2];
      Point p1 = g_vert[n1];
      Point p2 = g_vert[n2];
      Point p3 = g_vert[n3];
      Vector v1 = Vector(p1,p2);
      Vector v2 = Vector(p1,p3);
      Vector norm = (v1 % v2)*0.5;
      bx += norm[0];
      by += norm[1];
      bz += norm[2];
      wet[i] += sqrt(norm*norm);
    }
    fprintf(out_f,"\n%10d %14.7g %14.7g %14.7g %14.7g",i,bx,by,bz,wet[i]);
    if (i > 0)
    {
      tx += bx;
      ty += by;
      tz += bz;
      twet += wet[i];
    }
  }
  fprintf(out_f,"\n   Total   %14.7g %14.7g %14.7g %14.7g\n",tx,ty,tz,twet);
  fflush(out_f);

  // compute perimeter distance for each body, if exists
  double ds;
  for (c=0; c < n_chain; c++)
  {
    ds=0.0;
    for (i=1; i < chain[c].nn; i++)
    {
      n0 = chain[c].nodes[i-1];
      n1 = chain[c].nodes[i];
      ds += distance(g_vert[n0],g_vert[n1]);
    }
    g = chain[c].mat;
    perimeter[g] += ds;
  }

  Point lo = min_point();
  Point hi = max_point();
  layers = new int[ngb+1];
  g_space = new double[ngb+1];
  n_space = new double[ngb+1];
  for (i=0; i <= ngb; i++)
  {
    layers[i] = 0;
    double dx = bhi[i][0]-blo[i][0];
    double dy = bhi[i][1]-blo[i][1];
    double dz = bhi[i][2]-blo[i][2];
    double dmn = MIN(dx,MIN(dy,dz));
    double dmx = MAX(dx,MAX(dy,dz));
    double dmd = dx + dy + dz - dmn - dmx;
    // compute minimum of hydraulic diameter
    if (perimeter[i] > 1.0e-20)
      dmd = 4.0*wet[i]/perimeter[i];
    dmd *= 0.025;
    if (n_end[i]-n_begin[i] > 0)
    {
      dmx = MAX(hi[0]-lo[0],MAX(hi[1]-lo[1],hi[2]-lo[2]));
      while (dmx > dmd)
        dmx *= 0.5;
    } else
      dmx = dmd;
    g_space[i] = dmx;
    n_space[i] = dmx;
    for (j=n_begin[i]; j < n_end[i]; j++)
      g_facet[j].space = g_space[i];
  }
  delete[] wet;
  delete[] perimeter;

  // Check for existing spacing and layers parameter file
  sprintf(fname,"P_HUGG_geometry.params");
  if ((g_io=fopen(fname,"r")) == 0)
  {
    if (my_rank == 0)
    {
      // Create file with current values
      if ((g_io=fopen(fname,"w")) == 0)
      {
        fprintf(stderr,"\nCouldn't open parameter file for writing.");
        fflush(stderr);
        //MPI_Abort(MPI_COMM_WORLD,0);
        exit(0);
      }
      fprintf(g_io,"Number of boundary files\n");
      fprintf(g_io,"%d\n",ngf);
      for (i=0; i < ngf; i++)
        fprintf(g_io,"%s\n",gnames[i]);
      fprintf(g_io,"Number of boundaries\n");
      fprintf(g_io,"%d\n",ngb);
      fprintf(g_io,"Boundary #   Layers      g_space       n_space              boundary name\n");
      fprintf(g_io,"---------- ---------- -------------- -------------- --------------------------------\n");
      for (i=0; i <= ngb; i++)
      {
        if (i==0)
          sprintf(g_bname[i],"Virtual");
        else
          sprintf(g_bname[i],"Boundary_%d",i);
        fprintf(g_io,"%10d %10d %14.7g %14.7g %-32s\n",
                     i,layers[i],g_space[i],n_space[i],g_bname[i]);
      }
      fclose(g_io);
    }
  } else {
    int flag = 1;
    fgets(buff,bdim,g_io);
    fgets(buff,bdim,g_io);
    sscanf(buff,"%d",&j);
    if (j != ngf)
      flag = 0;
    for (i=0; i < j && flag; i++)
    {
      fgets(buff,bdim,g_io);
      if (strncmp(buff,gnames[i],strlen(gnames[i])) != 0)
        flag=0;
    }
    if (flag)
    {
      fgets(buff,bdim,g_io);
      fgets(buff,bdim,g_io);
      sscanf(buff,"%d",&i);
      if (i != ngb)
        flag = 0;
    }
    if (!flag)
    {
      if (my_rank == 0)
      {
        fprintf(out_f,"\nExisting parameter file does not contain the proper");
        fprintf(out_f,"\n number of boundaries.");
        fprintf(out_f,"\n Using default values for parameters!");
        fclose(g_io);
        // Create file with current values
        if ((g_io=fopen(fname,"w")) == 0)
        {
          fprintf(stderr,"\nCouldn't open parameter file for writing.");
          fflush(stderr);
          //MPI_Abort(MPI_COMM_WORLD,0);
          exit(0);
        }
        fprintf(g_io,"Number of boundary files\n");
        fprintf(g_io,"%d\n",ngf);
        for (i=0; i < ngf; i++)
          fprintf(g_io,"%s\n",gnames[i]);
        fprintf(g_io,"Number of boundaries\n");
        fprintf(g_io,"%d\n",ngb);
        fprintf(g_io,"Boundary #   Layers      g_space       n_space              boundary name\n");
        fprintf(g_io,"---------- ---------- -------------- -------------- --------------------------------\n");
        for (i=0; i <= ngb; i++)
        {
          if (i==0)
            sprintf(g_bname[i],"Virtual");
          else
            sprintf(g_bname[i],"Boundary_%d",i);
          fprintf(g_io,"%10d %10d %14.7g %14.7g %-32s\n",
                       i,layers[i],g_space[i],n_space[i],g_bname[i]);
        }
      }
    } else {
      fgets(buff,bdim,g_io);
      fgets(buff,bdim,g_io);
      for (i=0; i <= ngb; i++)
      {
        fgets(buff,bdim,g_io);
        sscanf(buff,"%d %d %lf %lf %32s",&j,&layers[i],&g_space[i],&n_space[i],g_bname[i]);
        //if (layers[i] > 0)
        //{
        //  double fac = 1.0;
        //  double ds = 1.15;
        //  for (j=0; j < 30; j++)
        //  {
        //    fac += ds;
        //    ds *= 1.15;
        //  }
        //  g_space[i] = MIN(g_space[i],n_space[i]*fac);
        //}
        for (j=n_begin[i]; j < n_end[i]; j++)
        {
          //if (layers[i])
          //  g_facet[j].space = MAX(g_facet[j].space,g_space[i]);
          //else
          //  g_facet[j].space = g_space[i];
          g_facet[j].space = g_space[i];
          g_facet[j].n_space = n_space[i];
        }
      }
    }
    fclose(g_io);
  }
  fprintf(out_f,"\nNumber of boundary files");
  fprintf(out_f,"\n%d",ngf);
  for (i=0; i < ngf; i++)
    fprintf(out_f,"\n%s",gnames[i]);
  fprintf(out_f,"\nNumber of boundaries");
  fprintf(out_f,"\n%d",ngb);
  fprintf(out_f,"\nBoundary #   Layers      g_space       n_space             boundary name");
  fprintf(out_f,"\n---------- ---------- -------------- -------------- --------------------------------");
  for (i=0; i <= ngb; i++)
  {
    fprintf(out_f,"\n%10d %10d %14.7g %14.7g %-32s",
                 i,layers[i],g_space[i],n_space[i],g_bname[i]);
  }
  fflush(out_f);

  int tenth;

  // modify facet spacing based on curvature with neighbors
  int (*nbr)[3];
  nbr = new int[n_gfacets][3];
  if (!nbr)
  {
    fprintf(out_f,"\nREAD_GEOM: unable to allocate memory for nbr!");
    fflush(out_f);
    //MPI_Abort(MPI_COMM_WORLD,0);
    exit(0);
  }
  for (i=0; i < n_gfacets; i++)
    nbr[i][0] = nbr[i][1] = nbr[i][2] = -1;
  Linked_List **Lhash;
  Lhash = (Linked_List**)malloc((n_gnodes)*sizeof(Linked_List*));
  if (!Lhash)
  {
    fprintf(out_f,"\nREAD_GEOM: unable to allocate memory for Lhash!");
    fflush(out_f);
    //MPI_Abort(MPI_COMM_WORLD,0);
    exit(0);
  }
  for (i=0; i < n_gnodes; i++)
  {
    Lhash[i] = new Linked_List();
    if (!Lhash[i])
    {
      fprintf(out_f,"\nREAD_GEOM: unable to allocate memory for Lhash[%d]!",i);
      fflush(out_f);
      //MPI_Abort(MPI_COMM_WORLD,0);
      exit(0);
    }
  }
  for (i=0; i < n_gfacets; i++)
  {
    n0 = g_facet[i].nodes[0];
    n1 = g_facet[i].nodes[1];
    n2 = g_facet[i].nodes[2];
    Lhash[n0]->Insert(i);
    Lhash[n1]->Insert(i);
    Lhash[n2]->Insert(i);
  }

  Linked_Node *hd0, *hd1;
  for (j=0; j < 3; j++)
  {
    for (i=0; i < n_gfacets; i++)
    {
      n0 = g_facet[i].nodes[j];
      if (j < 2)
        n1 = g_facet[i].nodes[j+1];
      else
        n1 = g_facet[i].nodes[0];
      hd0 = Lhash[n0]->head;
      m = -1;
      while (hd0 && m < 0)
      {
        n = hd0->data;
        if (n != i)
        {
          hd1 = Lhash[n1]->head;
          while (hd1 && m < 0)
          {
            if (hd1->data == n)
              m = n;
            hd1 = hd1->next;
          }
        }
        hd0 = hd0->next;
      }
      nbr[i][j] = m;
    }
  }
  for (i=0; i < n_gnodes; i++)
    delete Lhash[i];
  free(Lhash);

  int m0, m1, m2;
  double rad, lspace;
  Point cc;
  for (g=1; g <= ngb; g++)
  {
    if (layers[g] == 0 || layers[g] == 2)
      continue;
    tenth = (n_end[g]-n_begin[g]+1)/10;
    tenth = MAX(1,tenth);
    fprintf(out_f,"\nPerforming curvature test for boundary %d. Dot every %d facets.\n",g,tenth);
    fflush(out_f);
    for (i=n_begin[g]; i < n_end[g]; i++)
    {
      if ((i-n_begin[g]+1) % tenth == 0 || i == n_end[g]-1)
      {
        fprintf(out_f,".");
        fflush(out_f);
      }
      lspace = g_facet[i].space;
      n0 = g_facet[i].nodes[0];
      n1 = g_facet[i].nodes[1];
      n2 = g_facet[i].nodes[2];
      v1 = Vector(g_vert[n0],g_vert[n1]);
      v2 = Vector(g_vert[n0],g_vert[n2]);
      Vector norm = v1 % v2;
      norm.normalize();
      for (j=0; j < 3; j++)
      {
        switch(j)
        {
          case 0:
            n0 = g_facet[i].nodes[0];
            n1 = g_facet[i].nodes[1];
            n2 = g_facet[i].nodes[2];
            break;
          case 1:
            n0 = g_facet[i].nodes[1];
            n1 = g_facet[i].nodes[2];
            n2 = g_facet[i].nodes[0];
            break;
          case 2:
            n0 = g_facet[i].nodes[2];
            n1 = g_facet[i].nodes[0];
            n2 = g_facet[i].nodes[1];
            break;
        }
        if (nbr[i][j] < 0)
          continue;

        m0 = g_facet[nbr[i][j]].nodes[0];
        m1 = g_facet[nbr[i][j]].nodes[1];
        m2 = g_facet[nbr[i][j]].nodes[2];
        v1 = Vector(g_vert[m0],g_vert[m1]);
        v2 = Vector(g_vert[m0],g_vert[m2]);
        Vector v3 = v1 % v2;
        v3.normalize();
        
        rad = g_facet[i].space * pow((1.0+v3*norm)*0.5,2.0);
        //rad = circumcenter(g_vert[m0],g_vert[m1],g_vert[m2],g_vert[n2],cc);

        lspace = MIN(lspace,rad);
      }
      g_facet[i].space = MIN(g_facet[i].space,lspace);
    }
    fprintf(out_f,"\n");
    fflush(out_f);
  }
  delete[] nbr;

  // create Octree sorted facet lists
  Point p0, p1, p2;
  double olo[3], ohi[3];
  p0 = min_point();
  p1 = max_point();
  for (i=0; i < 3; i++)
  {
    olo[i] = p0[i];
    ohi[i] = p1[i];
  }

  Octree_Storage *dummy=0;
  Octree_root = new Octree_Storage((Octree_Storage*)dummy,olo,ohi);
  void* *ptr;
  double (*tlo)[3], (*thi)[3];
  ptr = new void*[n_gfacets];
  tlo = new double[n_gfacets][3];
  thi = new double[n_gfacets][3];
  for (i=0; i < n_gfacets; i++)
  {
    n0 = g_facet[i].nodes[0];
    n1 = g_facet[i].nodes[1];
    n2 = g_facet[i].nodes[2];
    tlo[i][0] = MIN(g_vert[n0][0],MIN(g_vert[n1][0],g_vert[n2][0]));
    tlo[i][1] = MIN(g_vert[n0][1],MIN(g_vert[n1][1],g_vert[n2][1]));
    tlo[i][2] = MIN(g_vert[n0][2],MIN(g_vert[n1][2],g_vert[n2][2]));
    thi[i][0] = MAX(g_vert[n0][0],MAX(g_vert[n1][0],g_vert[n2][0]));
    thi[i][1] = MAX(g_vert[n0][1],MAX(g_vert[n1][1],g_vert[n2][1]));
    thi[i][2] = MAX(g_vert[n0][2],MAX(g_vert[n1][2],g_vert[n2][2]));
    ptr[i] = (void*)&g_facet[i];
  }
  Octree_root->Store_In_Octree(n_gfacets,ptr,tlo,thi);
  delete[] ptr;
  delete[] tlo;
  delete[] thi;

  glist = new PList();
  glist->Redimension(n_gfacets);
  Octree_max_level = Octree_root->max_level(1);
  fprintf(out_f,"\nMaximum number of levels in Octree-sorted geometry= %d\n",Octree_max_level);
  fflush(out_f);

  double smn, smx;
  // perform intersection test 
  int gs, ge, dg;
  //int mdim;
  //int sposition, rposition, nreq_s, nreq_r, p;
  //char *sbuff, **rbuff;

  // allocate memory for MPI non-blocked info
  //int *recvcnt = new int[num_procs];
  //MPI_Request *srequest = new MPI_Request[num_procs];
  //MPI_Request *rrequest = new MPI_Request[num_procs];
  //MPI_Status *statuses  = new MPI_Status[num_procs];
  //mdim = (sizeof(double)+sizeof(int))*(n_gfacets/MAX(1,num_procs-1) + num_procs);
  //rbuff = new char*[num_procs];
  /*for (p=0; p < num_procs; p++)
  {
    rbuff[p] = new char[mdim];
  }
  sbuff = new char[mdim];*/

  for (g=1; g <= ngb; g++)
  {
    //MPI_Barrier(MPI_COMM_WORLD);

    if (layers[g] == 0 || layers[g] == 1)
      continue;
    fprintf(out_f,"\nPerforming intersection test for boundary %d.",g);

    dg = (n_end[g]-n_begin[g])/num_procs;
    gs = n_begin[g] + my_rank*dg;
    ge = gs + dg;
    if (my_rank == num_procs-1)
      ge = n_end[g];
    tenth = (ge-gs)/10;
    tenth = MAX(1,tenth);
    fprintf(out_f,"\nPerforming intersection/thinness test for geometry facets %d through %d. Dot every %d facets.\n",gs,ge-1,tenth);
    fflush(out_f);
    smn = 1.0e20;
    smx = -1.0e20;
    for (i=gs; i < ge; i++)
    {
      if ((i-gs+1) % tenth == 0 || i == ge-1)
      {
        fprintf(out_f,".");
        fflush(out_f);
      }
      n0 = g_facet[i].nodes[0];
      n1 = g_facet[i].nodes[1];
      n2 = g_facet[i].nodes[2];
      v1 = Vector(g_vert[n0],g_vert[n1]);
      v2 = Vector(g_vert[n0],g_vert[n2]);
      norm = v1 % v2;
      norm.normalize();
      double nudge = MIN(v1.magnitude(),v2.magnitude());
      v1 = Vector(g_vert[n1],g_vert[n2]);
      nudge = MIN(nudge,v1.magnitude())*0.01;
      p0 = (g_vert[n0]+g_vert[n1]+g_vert[n2])/3.0 + Point(norm[0],norm[1],norm[2])*nudge;
      p1 = p0 + Point(norm[0],norm[1],norm[2])*g_facet[i].n_space;
      if (Intersect(p0,p1,p2))
        g_facet[i].n_space = MIN(g_facet[i].n_space,distance(p0,p2)/2.0);
      p0 = (g_vert[n0]+g_vert[n1]+g_vert[n2])/3.0 - Point(norm[0],norm[1],norm[2])*nudge;
      p1 = p0 - Point(norm[0],norm[1],norm[2])*g_facet[i].n_space;
      if (Intersect(p0,p1,p2))
        g_facet[i].n_space = MIN(g_facet[i].n_space,distance(p0,p2)/2.0);
      smn = MIN(smn,g_facet[i].n_space);
      smx = MAX(smx,g_facet[i].n_space);
    }
    fprintf(out_f,"\nProcess %d minimum computed geometric normal spacing = %lg",my_rank,smn);
    fprintf(out_f,"\nProcess %d maximum computed geometric normal spacing = %lg",my_rank,smx);
    fprintf(out_f,"\n");
    fflush(out_f);

    // Exchange spacing values with other processors
    //MPI_Barrier(MPI_COMM_WORLD);

    /*int sendcnt = ge-gs;

    nreq_s = nreq_r = 0;
    for (p=0; p < num_procs; p++)
    {
      recvcnt[p] = sendcnt;
      if (p != my_rank)
      {
        MPI_Isend(&sendcnt,1,MPI_INT,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
        MPI_Irecv(&recvcnt[p],1,MPI_INT,p,      p,MPI_COMM_WORLD,&rrequest[nreq_r]);
        nreq_s++;
        nreq_r++;
      }
    }
    //MPI_Waitall(nreq_s,srequest,statuses);
    //MPI_Waitall(nreq_r,rrequest,statuses);

    sposition = 0;
    for (i=gs; i < ge; i++)
    {
      MPI_Pack(&i,1,MPI_INT,sbuff,mdim,&sposition,MPI_COMM_WORLD);
      MPI_Pack(&g_facet[i].n_space,1,MPI_DOUBLE,sbuff,mdim,&sposition,MPI_COMM_WORLD);
    }

    // exchange messages
    nreq_s = nreq_r = 0;
    for (p=0; p < num_procs; p++)
    {
      if (p != my_rank)
      {
        MPI_Isend(sbuff,mdim,MPI_CHAR,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
        MPI_Irecv(rbuff[p],mdim,MPI_CHAR,p,      p,MPI_COMM_WORLD,&rrequest[nreq_r]);
        nreq_s++;
        nreq_r++;
      }
    }
    MPI_Waitall(nreq_s,srequest,statuses);
    MPI_Waitall(nreq_r,rrequest,statuses);

    for (p=0; p < num_procs; p++)
    {
      if (p == my_rank)
        continue;

      rposition=0;
      for (i=0; i < recvcnt[p]; i++)
      {
        MPI_Unpack(rbuff[p],mdim,&rposition,&j,1,MPI_INT,MPI_COMM_WORLD);
        MPI_Unpack(rbuff[p],mdim,&rposition,&ds,1,MPI_DOUBLE,MPI_COMM_WORLD);
        g_facet[j].n_space = ds;
      }
    }*/
  }

  /*for (p=0; p < num_procs; p++)
    delete[] rbuff[p];
  delete[] rbuff;
  delete[] sbuff;
  delete[] statuses;
  delete[] rrequest;
  delete[] srequest;
  delete[] recvcnt;*/

  //double *stmp = new double[n_gfacets];
  //for (j=0; j < 3; j++)
  //{
  //  for (i=0; i < n_gfacets; i++)
  //    stmp[i] = g_facet[i].space;
  //  for (i=0; i < n_gfacets; i++)
  //  {
  //    n0 = nbr[i][0];
  //    if (n0 < 0) n0 = i;
  //    n1 = nbr[i][1];
  //    if (n1 < 0) n1 = i;
  //    n2 = nbr[i][2];
  //    if (n2 < 0) n2 = i;
  //    ds = (stmp[n0]+stmp[n1]+stmp[n2])/3.0;
  //    ds = MIN(0.0,ds-g_facet[i].space);
  //    g_facet[i].space += ds*0.25;
  //  }
  //}
  //delete[] stmp;

  for (g=1; g <= ngb; g++)
  {
    fprintf(out_f,"\nOverall spacing limits for boundary %d:",g);
    smn = 1.0e20;
    smx = -1.0e20;
    for (i=n_begin[g]; i < n_end[g]; i++)
    {
      smn = MIN(smn,g_facet[i].space);
      smx = MAX(smx,g_facet[i].space);
    }
    fprintf(out_f,"\n	Minimum geometric tangential spacing = %lg",smn);
    fprintf(out_f,"\n	Maximum geometric tangential spacing = %lg",smx);
    smn = 1.0e20;
    smx = -1.0e20;
    for (i=n_begin[g]; i < n_end[g]; i++)
    {
      smn = MIN(smn,g_facet[i].n_space);
      smx = MAX(smx,g_facet[i].n_space);
    }
    fprintf(out_f,"\n	Minimum geometric normal spacing = %lg",smn);
    fprintf(out_f,"\n	Maximum geometric normal spacing = %lg",smx);
  }
  fflush(out_f);

#ifdef _DEBUG
  Plotfile(1);
#endif

}

int geometry::Identify_material_combos()
{
  int i, j, k, g, m, n;
  int n0, n1, n2;

  // Initialize with list of single material states
  n_combos = 0;
  mat_combo = (Mat_Combo*)malloc(ngb*sizeof(Mat_Combo));
  for (i=0; i < ngb; i++)
  {
    mat_combo[i].nm = 0;
    mat_combo[i].mat = 0;
    mat_combo[i].Add_mat(i+1);
  }
  n_combos = ngb;
  
  int **map;
  map = new int*[n_gnodes];
  for (n=0; n < n_gnodes; n++)
  {
    map[n] = new int[ngb];
    for (j=0; j < ngb; j++)
      map[n][j] = 0;
  }

  // generate table for each geometry node
  // that indicates which material it is associated with
  for (g=1; g <= ngb; g++)
  {
    for (i=n_begin[g]; i < n_end[g]; i++)
    {
      n0 = g_facet[i].nodes[0];
      n1 = g_facet[i].nodes[1];
      n2 = g_facet[i].nodes[2];
      map[n0][g-1] = 1;
      map[n1][g-1] = 1;
      map[n2][g-1] = 1;
    }
  }

  int mflag,iflag;
  // check each node for new material combinations
  for (n=0; n < n_gnodes; n++)
  {
    for (k=i=0; i < ngb; i++)
      k += map[n][i];
    if (k < 2)
      continue;
    // compare against current states
    mflag = 0;
    for (m=ngb; m < n_combos && !mflag; m++)
    {
      iflag = 1;
      for (i=0; i < mat_combo[m].nm && iflag; i++)
      {
        j = mat_combo[m].mat[i]-1;
        if (map[n][j] == 0)
          iflag = 0;
      }
      mflag = iflag;
    }
    if (!mflag)
    {
      mat_combo=(Mat_Combo*)realloc((void*)mat_combo,(n_combos+1)*sizeof(Mat_Combo));
      mat_combo[n_combos].nm = 0;
      mat_combo[n_combos].mat = 0;
      for (i=0; i < ngb; i++)
        if (map[n][i] == 1)
          mat_combo[n_combos].Add_mat(i+1);
      n_combos++;
    }
  }

  // add other possible two material states
  k=n_combos;
  for (g=ngb; g < k; g++)
  {
    if (mat_combo[g].nm < 3)
      continue;
    for (i=0; i < mat_combo[g].nm-1; i++)
    {
      for (j=i+1; j < mat_combo[g].nm; j++)
      {
        mflag = 0;
        for (m=ngb; m < n_combos && !mflag; m++)
        {
          if (mat_combo[m].nm != 2)
            continue;
          if (mat_combo[g].mat[i] == mat_combo[m].mat[0] &&
              mat_combo[g].mat[j] == mat_combo[m].mat[1])
            mflag = 1;
        }
        if (!mflag)
        {
          mat_combo=(Mat_Combo*)realloc((void*)mat_combo,(n_combos+1)*sizeof(Mat_Combo));
          mat_combo[n_combos].nm = 0;
          mat_combo[n_combos].mat = 0;
          mat_combo[n_combos].Add_mat(mat_combo[g].mat[i]);
          mat_combo[n_combos].Add_mat(mat_combo[g].mat[j]);
          n_combos++;
        }
      }
    }
  }

  for (i=0; i < n_gnodes; i++)
    delete[] map[i];
  delete[] map;

  return(n_combos);

}

Point geometry::closest_chain(Point p, int &b)
{
  int c, i, n0, n1;
  Point pp, pc;
  double ds;
  double dist = 1.0e20;
  pc = p;
  for (c=0; c < n_chain; c++)
  {
    if (b > 0 && chain[c].mat != b)
      continue;
    for (i=1; i < chain[c].nn; i++)
    {
      n0 = chain[c].nodes[i-1];
      n1 = chain[c].nodes[i];
      ds=closest_to_edge(g_vert[n0],g_vert[n1],p,pp);
      if (ds < dist)
      {
        pc = pp;
        dist = ds;
      }
    }
  }
  return (pc);
}

Point geometry::closest(Point p, int &b, Vector &nrm)
{
  int g, i, l, n0, n1, n2, tb;
  Point p0, p1, p2, flo, fhi, ptemp, ps, ptol;
  Vector v1, v2, norm, tnrm;
  double mag, dist, dot;
  double pt[3], tol[3];
  Geom_Facet *gptr, *ptr;

  tb = b;
  dist = 1.0e20;
  gptr = NULL;
  b = 0;
  tnrm = nrm;
  tnrm.normalize();

  for (g=1; g <= ngb; g++)
  {
    if (tb > 0 && g != tb)
      continue;
    v1 = Vector(blo[g],bhi[g])*0.5;
    ptemp = (blo[g]+bhi[g])*0.5;
    v2 = Vector(p,ptemp);
    mag = v1.magnitude()+v2.magnitude();
    dist = MIN(dist,mag);
  }
  ptol = Point(dist,dist,dist);

  int mxlvl = Octree_max_level;

  pt[0] = p[0];
  pt[1] = p[1];
  pt[2] = p[2];
  for (l=1; l <= mxlvl; l++)
  {
    tol[0] = ptol[0];
    tol[1] = ptol[1];
    tol[2] = ptol[2];
    glist->max = 0;
    Octree_root->retrieve_list(pt,tol,l,glist);

    for (i=0; i < glist->max; i++)
    {
      ptr=(Geom_Facet*)glist->list[i];

      if (ptr->mat <= 0 || (tb > 0 && ptr->mat != tb))
        continue;
      n0 = ptr->nodes[0];
      n1 = ptr->nodes[1];
      n2 = ptr->nodes[2];
      p0 = g_vert[n0];
      p1 = g_vert[n1];
      p2 = g_vert[n2];

      // perform extent test using current closest distance
      flo[0] = MIN(p0[0],MIN(p1[0],p2[0]))-dist;
      flo[1] = MIN(p0[1],MIN(p1[1],p2[1]))-dist;
      flo[2] = MIN(p0[2],MIN(p1[2],p2[2]))-dist;
      fhi[0] = MAX(p0[0],MAX(p1[0],p2[0]))+dist;
      fhi[1] = MAX(p0[1],MAX(p1[1],p2[1]))+dist;
      fhi[2] = MAX(p0[2],MAX(p1[2],p2[2]))+dist;
      if (p[0] < flo[0] || p[0] > fhi[0] ||
          p[1] < flo[1] || p[1] > fhi[1] ||
          p[2] < flo[2] || p[2] > fhi[2])
        continue;

      mag = closest_to_facet(p0,p1,p2,p,ptemp);

      v1 = Vector(p0,p1);
      v2 = Vector(p0,p2);
      norm = v1 % v2;
      norm.normalize();
      dot = tnrm * norm;

      if (mag < dist && dot > -0.001)
      {
        dist = mag;
        ps = ptemp;
        gptr = ptr;
        nrm = norm;
        ptol = Point(dist*1.01,dist*1.01,dist*1.01);
      }
    }
  }

  if (tb <= 0 && gptr != NULL)
    b = gptr->mat;
  if (tb > 0)
    b = tb;

  return(ps);
}

Point geometry::min_point()
{
// Return minimum values of edge vertex coordinates
  Point p(1.0e20, 1.0e20, 1.0e20);
  for (int i=0; i < n_gnodes; i++)
  {
    if (g_vert[i][0] < p[0])
      p[0] = g_vert[i][0];
    if (g_vert[i][1] < p[1])
      p[1] = g_vert[i][1];
    if (g_vert[i][2] < p[2])
      p[2] = g_vert[i][2];
  }
  return p;
}

Point geometry::max_point()
{
// Return maximum values of edge vertex coordinates
  Point p(-1.0e20, -1.0e20, -1.0e20);
  for (int i=0; i < n_gnodes; i++)
  {
    if (g_vert[i][0] > p[0])
      p[0] = g_vert[i][0];
    if (g_vert[i][1] > p[1])
      p[1] = g_vert[i][1];
    if (g_vert[i][2] > p[2])
      p[2] = g_vert[i][2];
  }
  return p;
}

int geometry::Intersect(Point pa, Point pb, Point &p)
{
  int flag, i, l, n0, n1, n2;
  double ds, dmin;
  Point glo, ghi, lo, mid, hi, p0, p1, p2, ptemp, ptol;
  double pt[3],tol[3];
  Geom_Facet *ptr;

  lo[0] = MIN(pa[0],pb[0]);
  lo[1] = MIN(pa[1],pb[1]);
  lo[2] = MIN(pa[2],pb[2]);
  hi[0] = MAX(pa[0],pb[0]);
  hi[1] = MAX(pa[1],pb[1]);
  hi[2] = MAX(pa[2],pb[2]);
  mid = (lo + hi)*0.5;
  double eps = 1.0e-15;

  ptol = Point(MAX(hi[0]-lo[0],eps),
               MAX(hi[1]-lo[1],eps),
               MAX(hi[2]-lo[2],eps))*0.51;

  pt[0] = mid[0];
  pt[1] = mid[1];
  pt[2] = mid[2];
  tol[0] = ptol[0];
  tol[1] = ptol[1];
  tol[2] = ptol[2];
  glist->max = 0;
  l = -1;
  Octree_root->retrieve_list(pt,tol,l,glist);

  flag = 0;
  dmin = 1.0e20;
  for (i=0; i < glist->max; i++)
  {
    ptr = (Geom_Facet*)glist->list[i];
    if (ptr->mat <= 0)
      continue;
    n0 = ptr->nodes[0];
    n1 = ptr->nodes[1];
    n2 = ptr->nodes[2];
    p0 = g_vert[n0];
    p1 = g_vert[n1];
    p2 = g_vert[n2];
    glo[0] = MIN(p0[0],MIN(p1[0],p2[0]));
    glo[1] = MIN(p0[1],MIN(p1[1],p2[1]));
    glo[2] = MIN(p0[2],MIN(p1[2],p2[2]));
    ghi[0] = MAX(p0[0],MAX(p1[0],p2[0]));
    ghi[1] = MAX(p0[1],MAX(p1[1],p2[1]));
    ghi[2] = MAX(p0[2],MAX(p1[2],p2[2]));
    if (lo[0]-eps > ghi[0] || lo[1]-eps > ghi[1] || lo[2]-eps > ghi[2] ||
        hi[0]+eps < glo[0] || hi[1]+eps < glo[1] || hi[2]+eps < glo[2])
      continue;

    if (line_facet_intersect(pa, pb, p0, p1, p2, ptemp, eps))
    {
      flag = 1;
      ds = distance(pa,ptemp);
      if (ds < dmin)
      {
        dmin = ds;
        p = ptemp;
        lo[0] = MIN(pa[0],p[0]);
        lo[1] = MIN(pa[1],p[1]);
        lo[2] = MIN(pa[2],p[2]);
        hi[0] = MAX(pa[0],p[0]);
        hi[1] = MAX(pa[1],p[1]);
        hi[2] = MAX(pa[2],p[2]);
      }
    }
  }

  return(flag);

}

// Fieldview routines
#include "FV_routines.h"

void geometry::Plotfile(int mode)
{
  FILE *fp;
  char filename[80], buff[80];
  int g, i, n, csize, isize, fsize;
  int ibuf[7], *bnum;
  float *f;

  // Write geometry out as Fieldview file

  fp = 0;
  csize = sizeof(char);
  isize = sizeof(int);
  fsize = sizeof(float);

  bnum = new int[ngb+1];

  // create object name entry for later use
  sprintf(filename,"P_HUGG_geom_%d",my_rank);

  switch(mode)
  {
    case 0: // Binary
      strcat(filename,".uns");
      // Open file for write
      if ((fp = fopen(filename,"w")) == 0)
      {
        fprintf(out_f,"\nError opening file <%s>.",filename);
        fflush(out_f);
        //MPI_Abort(MPI_COMM_WORLD,0);
        exit(0);
      }

      // Fieldview bit pattern
      ibuf[0] = FV_MAGIC;
      fwrite(ibuf,isize,1,fp);

      // name and version number
      sprintf(buff,"FIELDVIEW");
      fwrite(buff,csize,80,fp);
      ibuf[0] = 3;
      ibuf[1] = 0;
      fwrite(ibuf,isize,2,fp);
      ibuf[0] = FV_GRIDS_FILE;
      fwrite(ibuf,isize,1,fp);
      ibuf[0] = 0;
      fwrite(ibuf,isize,1,fp);

      // Number of grids
      ibuf[0] = 1;
      fwrite(ibuf,isize,1,fp);

      // Number of boundary types
      ibuf[0] = ngb;
      fwrite(ibuf,isize,1,fp);

      // Boundary types
      for (g=1; g <= ngb; g++)
      {
        sprintf(buff,"Boundary %d",g);
        ibuf[0]=0;
        ibuf[1]=1;
        fwrite(ibuf,isize,2,fp);
        fwrite(ibuf,csize,80,fp);
      }

      // Node information  
      // (writing double points to create flat prism volume cells)
      ibuf[0] = FV_NODES;
      ibuf[1] = n_gnodes*2;
      fwrite(ibuf,isize,2,fp);
      //f = (float*)malloc(n_gnodes*2*fsize);
      f = new float[n_gnodes*2];
      // X coordinate
      for (i=n=0; n < n_gnodes; n++)
        f[i++] = (float)g_vert[n][0];
      for (n=0; n < n_gnodes; n++)
        f[i++] = (float)g_vert[n][0];
      fwrite(f,fsize,n_gnodes*2,fp);
      // Y coordinate
      for (i=n=0; n < n_gnodes; n++)
        f[i++] = (float)g_vert[n][1];
      for (n=0; n < n_gnodes; n++)
        f[i++] = (float)g_vert[n][1];
      fwrite(f,fsize,n_gnodes*2,fp);
      // Z coordinate
      for (i=n=0; n < n_gnodes; n++)
        f[i++] = (float)g_vert[n][2];
      for (n=0; n < n_gnodes; n++)
        f[i++] = (float)g_vert[n][2];
      fwrite(f,fsize,n_gnodes*2,fp);

      //free(f);
      delete[] f;

      // Boundary section
      for (g=1; g <= ngb; g++)
      {
        ibuf[0] = FV_FACES;
        ibuf[1] = g;
        ibuf[2] = n_end[g]-n_begin[g]+1;
        fwrite(ibuf,isize,3,fp);
        ibuf[3] = 0;
        for (int c=n_begin[g]; c < n_end[g]; c++)
        {
          ibuf[0] = g_facet[c].nodes[0]+1;
          ibuf[1] = g_facet[c].nodes[1]+1;
          ibuf[2] = g_facet[c].nodes[2]+1;
          fwrite(ibuf,isize,4,fp);
        }
      }

      // Element section
      int walls[5];
      walls[0] = NOT_A_WALL;
      walls[1] = NOT_A_WALL;
      walls[2] = NOT_A_WALL;
      walls[3] = A_WALL;
      walls[4] = NOT_A_WALL;
      ibuf[0] = FV_ELEMENTS;
      ibuf[1] = 0;          // # of tetrahedron
      ibuf[2] = 0;          // # of hexahedron
      ibuf[3] = n_gfacets;  // # of prism
      ibuf[4] = 0;          // # of pyramid
      fwrite(ibuf,isize,5,fp);
      for (g=0; g < n_gfacets; g++)
      {
        ibuf[0]=fv_encode_elem_header(FV_PRISM_ELEM_ID,walls);
        ibuf[1]=g_facet[g].nodes[0]+1;
        ibuf[2]=g_facet[g].nodes[0]+1+n_gnodes;
        ibuf[3]=g_facet[g].nodes[1]+1+n_gnodes;
        ibuf[4]=g_facet[g].nodes[1]+1;
        ibuf[5]=g_facet[g].nodes[2]+1+n_gnodes;
        ibuf[6]=g_facet[g].nodes[2]+1;
        fwrite(ibuf,isize,7,fp);
      }

      break;
    case 1:  // ASCII
      strcat(filename,"_ascii.uns");
      fprintf(out_f,"\nFilename = <%s>",filename);
      // Open file for write
      if ((fp = fopen(filename,"w")) == 0)
      {
        fprintf(out_f,"\nError opening file <%s>.",filename);
        fflush(out_f);
        //MPI_Abort(MPI_COMM_WORLD,0);
        exit(0);
      }
      // name and version number
      fprintf(fp,"FIELDVIEW 2 4\n");

      fprintf(fp,"Constants\n");
      //f = (float*)malloc(4*fsize);
      f = new float[4];
      // Time
      f[0] = 0.0;
      // Mach
      f[1] = 0.0;
      // Alpha
      f[2] = 0.0;
      // Reynold's number
      f[3] = 0.0;
      fprintf(fp,"%f %f %f %f\n",f[0],f[1],f[2],f[3]);
      //free(f);
      delete[] f;

      // Number of grids
      fprintf(fp,"Grids 1\n");

      // Number of boundary types
      i=0;
      for (n=0; n <= ngb; n++)
      {
        bnum[n]=0;
        for (g=0; g < n_gfacets; g++)
          if (g_facet[g].mat == n)
            bnum[n]++;
        if (bnum[n])
          i++;
      }
      fprintf(fp,"Boundary Table %d\n",i);

      // Boundary types
      for (g=1; g <= ngb; g++)
      {
        if (bnum[g])
          fprintf(fp,"1 0 1 Boundary %d\n",g);
      }
      if (bnum[0])
        fprintf(fp,"1 0 1 Boundary %d\n",ngb+1);

      // Number of variables
      fprintf(fp,"Variable Names 0\n");

      // Node information  
      // (writing double points to create flat prism volume cells)
      fprintf(fp,"Nodes %d\n",n_gnodes*2);
      for (n=0; n < n_gnodes; n++)
        fprintf(fp,"%lf %lf %lf\n",g_vert[n][0],g_vert[n][1],g_vert[n][2]);
      for (n=0; n < n_gnodes; n++)
        fprintf(fp,"%lf %lf %lf\n",g_vert[n][0],g_vert[n][1],g_vert[n][2]);

      // Boundary section
      fprintf(fp,"Boundary Faces %d\n",n_gfacets);
      for (g=0; g <= ngb; g++)
      {
        if (bnum[g])
        {
          for (int c=0; c < n_gfacets; c++)
          {
            if (g_facet[c].mat == g)
            {
              ibuf[0] = g_facet[c].nodes[0]+1;
              ibuf[1] = g_facet[c].nodes[1]+1;
              ibuf[2] = g_facet[c].nodes[2]+1;
              if (g == 0)
                fprintf(fp,"%d 3 %d %d %d\n",ngb+1,ibuf[0],ibuf[1],ibuf[2]);
              else
                fprintf(fp,"%d 3 %d %d %d\n",g,ibuf[0],ibuf[1],ibuf[2]);
            }
          }
        }
      }

      // Element section
      fprintf(fp,"Elements\n");
      for (g=0; g < n_gfacets; g++)
      {
        ibuf[0]=g_facet[g].nodes[0]+1;
        ibuf[1]=g_facet[g].nodes[0]+1+n_gnodes;
        ibuf[2]=g_facet[g].nodes[1]+1+n_gnodes;
        ibuf[3]=g_facet[g].nodes[1]+1;
        ibuf[4]=g_facet[g].nodes[2]+1+n_gnodes;
        ibuf[5]=g_facet[g].nodes[2]+1;
        fprintf(fp,"3 1\n %d %d %d %d %d %d\n",ibuf[0],ibuf[1],ibuf[2],
                                               ibuf[3],ibuf[4],ibuf[5]);
      }

      // Variable section
      fprintf(fp,"Variables\n");

      // Boundary variable section
      fprintf(fp,"Boundary Variables\n");
      break;
    default:
      break;
  }

  fclose(fp);

  delete[] bnum;
}

int geometry::In_Out(Point pt, double tol)
{
  int flag, n0, n1, n2;
  Vector dv, dp;
  Point *ept, *wpt, *npt, *spt, *upt, *dpt;
  Point po, pp, p0, p1, p2;
  int i, j, g, ebp, wbp, nbp, sbp, ubp, dbp;
  int edim, wdim, ndim, sdim, udim, ddim;

  // Use ray tracing to determine in/out status of test point

  // Initially dimension boundary intersection arrays by the number of
  // geometry bodies. Later this may need to be expanded.
  ept=(Point*)malloc(ngb*sizeof(Point));
  wpt=(Point*)malloc(ngb*sizeof(Point));
  npt=(Point*)malloc(ngb*sizeof(Point));
  spt=(Point*)malloc(ngb*sizeof(Point));
  upt=(Point*)malloc(ngb*sizeof(Point));
  dpt=(Point*)malloc(ngb*sizeof(Point));
  edim = wdim = ndim = sdim = udim = ddim = ngb;
  
  dp = Vector(min_point(),max_point())*10.0;

  // east ray trace
  ebp = 0;
  po = pt;
  po[0] += dp[0];

  for (g=0; g < n_gfacets; g++)
  {
    if (g_facet[g].mat <= 0)
      continue;
    n0 = g_facet[g].nodes[0];
    n1 = g_facet[g].nodes[1];
    n2 = g_facet[g].nodes[2];
    p0 = g_vert[n0];
    p1 = g_vert[n1];
    p2 = g_vert[n2];
    if (MIN(p0[1],MIN(p1[1],p2[1])) > pt[1] ||
        MAX(p0[1],MAX(p1[1],p2[1])) < pt[1] ||
        MIN(p0[2],MIN(p1[2],p2[2])) > pt[2] ||
        MAX(p0[2],MAX(p1[2],p2[2])) < pt[2])
      continue;
    if (line_facet_intersect(po,pt,p0,p1,p2,pp,tol))
    {
      j = -1;
      for (i=0; i < ebp && j == -1; i++)
      {
        if (distance(ept[i],pp) < tol)
          j = i;
      }
      if (j == -1)
      {
        if (ebp == edim)
        {
          edim += ngb;
          ept=(Point*)realloc((void*)ept,edim*sizeof(Point));
        }
        ept[ebp++] = pp;
      }
    }
  }
        
  // west ray trace
  wbp = 0;
  po = pt;
  po[0] -= dp[0];

  for (g=0; g < n_gfacets; g++)
  {
    if (g_facet[g].mat <= 0)
      continue;
    n0 = g_facet[g].nodes[0];
    n1 = g_facet[g].nodes[1];
    n2 = g_facet[g].nodes[2];
    p0 = g_vert[n0];
    p1 = g_vert[n1];
    p2 = g_vert[n2];
    if (MIN(p0[1],MIN(p1[1],p2[1])) > pt[1] ||
        MAX(p0[1],MAX(p1[1],p2[1])) < pt[1] ||
        MIN(p0[2],MIN(p1[2],p2[2])) > pt[2] ||
        MAX(p0[2],MAX(p1[2],p2[2])) < pt[2])
      continue;
    if (line_facet_intersect(po,pt,p0,p1,p2,pp,tol))
    {
      j = -1;
      for (i=0; i < wbp && j == -1; i++)
      {
        if (distance(wpt[i],pp) < tol)
          j = i;
      }
      if (j == -1)
      {
        if (wbp == wdim)
        {
          wdim += ngb;
          wpt=(Point*)realloc((void*)wpt,wdim*sizeof(Point));
        }
        wpt[wbp++] = pp;
      }
    }
  }
        
  // south ray trace
  sbp = 0;
  po = pt;
  po[1] -= dp[1];

  for (g=0; g < n_gfacets; g++)
  {
    if (g_facet[g].mat <= 0)
      continue;
    n0 = g_facet[g].nodes[0];
    n1 = g_facet[g].nodes[1];
    n2 = g_facet[g].nodes[2];
    p0 = g_vert[n0];
    p1 = g_vert[n1];
    p2 = g_vert[n2];
    if (MIN(p0[0],MIN(p1[0],p2[0])) > pt[0] ||
        MAX(p0[0],MAX(p1[0],p2[0])) < pt[0] ||
        MIN(p0[2],MIN(p1[2],p2[2])) > pt[2] ||
        MAX(p0[2],MAX(p1[2],p2[2])) < pt[2])
      continue;
    if (line_facet_intersect(po,pt,p0,p1,p2,pp,tol))
    {
      j = -1;
      for (i=0; i < sbp && j == -1; i++)
      {
        if (distance(spt[i],pp) < tol)
          j = i;
      }
      if (j == -1)
      {
        if (sbp == sdim)
        {
          sdim += ngb;
          spt=(Point*)realloc((void*)spt,sdim*sizeof(Point));
        }
        spt[sbp++] = pp;
      }
    }
  }
        
  // north ray trace
  nbp = 0;
  po = pt;
  po[1] += dp[1];

  for (g=0; g < n_gfacets; g++)
  {
    if (g_facet[g].mat <= 0)
      continue;
    n0 = g_facet[g].nodes[0];
    n1 = g_facet[g].nodes[1];
    n2 = g_facet[g].nodes[2];
    p0 = g_vert[n0];
    p1 = g_vert[n1];
    p2 = g_vert[n2];
    if (MIN(p0[0],MIN(p1[0],p2[0])) > pt[0] ||
        MAX(p0[0],MAX(p1[0],p2[0])) < pt[0] ||
        MIN(p0[2],MIN(p1[2],p2[2])) > pt[2] ||
        MAX(p0[2],MAX(p1[2],p2[2])) < pt[2])
      continue;
    if (line_facet_intersect(po,pt,p0,p1,p2,pp,tol))
    {
      j = -1;
      for (i=0; i < nbp && j == -1; i++)
      {
        if (distance(npt[i],pp) < tol)
          j = i;
      }
      if (j == -1)
      {
        if (nbp == ndim)
        {
          ndim += ngb;
          npt=(Point*)realloc((void*)npt,ndim*sizeof(Point));
        }
        npt[nbp++] = pp;
      }
    }
  }

  // down ray trace
  dbp = 0;
  po = pt;
  po[2] -= dp[2];

  for (g=0; g < n_gfacets; g++)
  {
    if (g_facet[g].mat <= 0)
      continue;
    n0 = g_facet[g].nodes[0];
    n1 = g_facet[g].nodes[1];
    n2 = g_facet[g].nodes[2];
    p0 = g_vert[n0];
    p1 = g_vert[n1];
    p2 = g_vert[n2];
    if (MIN(p0[0],MIN(p1[0],p2[0])) > pt[0] ||
        MAX(p0[0],MAX(p1[0],p2[0])) < pt[0] ||
        MIN(p0[1],MIN(p1[1],p2[1])) > pt[1] ||
        MAX(p0[1],MAX(p1[1],p2[1])) < pt[1])
      continue;
    if (line_facet_intersect(po,pt,p0,p1,p2,pp,tol))
    {
      j = -1;
      for (i=0; i < dbp && j == -1; i++)
      {
        if (distance(dpt[i],pp) < tol)
          j = i;
      }
      if (j == -1)
      {
        if (dbp == ddim)
        {
          ddim += ngb;
          dpt=(Point*)realloc((void*)dpt,ddim*sizeof(Point));
        }
        dpt[dbp++] = pp;
      }
    }
  }

  // up ray trace
  ubp = 0;
  po = pt;
  po[2] += dp[2];

  for (g=0; g < n_gfacets; g++)
  {
    if (g_facet[g].mat <= 0)
      continue;
    n0 = g_facet[g].nodes[0];
    n1 = g_facet[g].nodes[1];
    n2 = g_facet[g].nodes[2];
    p0 = g_vert[n0];
    p1 = g_vert[n1];
    p2 = g_vert[n2];
    if (MIN(p0[0],MIN(p1[0],p2[0])) > pt[0] ||
        MAX(p0[0],MAX(p1[0],p2[0])) < pt[0] ||
        MIN(p0[1],MIN(p1[1],p2[1])) > pt[1] ||
        MAX(p0[1],MAX(p1[1],p2[1])) < pt[1])
      continue;
    if (line_facet_intersect(po,pt,p0,p1,p2,pp,tol))
    {
      j = -1;
      for (i=0; i < ubp && j == -1; i++)
      {
        if (distance(upt[i],pp) < tol)
          j = i;
      }
      if (j == -1)
      {
        if (ubp == udim)
        {
          udim += ngb;
          upt=(Point*)realloc((void*)upt,udim*sizeof(Point));
        }
        upt[ubp++] = pp;
      }
    }
  }

  // now vote
  int iflag, oflag;
  iflag = oflag = 0;
  if ((ebp % 2) == 0)
    oflag++;
  else
    iflag++;
  if ((wbp % 2) == 0)
    oflag++;
  else
    iflag++;
  if ((sbp % 2) == 0)
    oflag++;
  else
    iflag++;
  if ((nbp % 2) == 0)
    oflag++;
  else
    iflag++;
  if ((dbp % 2) == 0)
    oflag++;
  else
    iflag++;
  if ((ubp % 2) == 0)
    oflag++;
  else
    iflag++;

  if (oflag == 6)
    flag = -1;
  else if (iflag == 6)
    flag = 1;
  else
    flag = 0;

  free(ept);
  free(wpt);
  free(npt);
  free(spt);
  free(dpt);
  free(upt);

  return(flag);

}
