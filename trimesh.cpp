#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/param.h>
#include <assert.h>
#include "Linked_List.h"
#include "Point2D.h"
#include "Vector2D.h"
#include "trimesh.h"

#define ABS(x) ((x) >= 0 ? (x) : -(x))
#define SIGN(a,b) ((b) < 0.0 ? -ABS(a) : ABS(a))

int circle_test(Point2D p1, Point2D p2, Point2D p3, Point2D p);
double max_angle(Point2D p1, Point2D p2, Point2D p3);
double dist( Point2D p1, Point2D p2);
void make_nbrs(Linked_List *mknbr, int nn, int nt, int tri[][3], int nbr[][3], Linked_List **Lhash);
void nbr_search(Point2D pt, Point2D p[], int tri[][3], int nbr[][3], Linked_List *del, int m);
void compress_list(Linked_List *del, Linked_List **Lhash, Linked_List *mknbr, int tri[][3], int nbr[][3], int &nt);
int search(int t, Point2D pt, Point2D p[], int tri[][3], int nbr[][3]);
double det_3X3(double m[3][3]);
//int trimesh(const int npt, const int tdim, const int nb, int nab, int *nbs, int ***bs, Point2D pt[], int tri[][3]);

// Wrapper function for use with individual coordinate arrays
int trimesh(const int npt, const int tdim, const int nb, int nab, int *nbs, int ***bs, double x[], double y[], int tri[][3])
{
  int nt;
  Point2D *pt;

  pt = new Point2D[npt];
  for (nt=0; nt < npt; nt++)
    pt[nt] = Point2D(x[nt],y[nt]);

  nt = trimesh(npt,tdim,nb,nab,nbs,bs,pt,tri);

  delete[] pt;

  return(nt);
}

int trimesh(const int npt, const int tdim, const int nb, int nab, int *nbs, int ***bs, Point2D pt[], int tri[][3])
{
  int b, i, j, m, n, nn, n0, n1, n2, t;
  int seed, m0, m1, m2;
  Linked_List **Lhash;
  Linked_List *mknbr;
  Linked_Node *hd;
  
  //for (i = 0; i < nbs[0]; i++)
    //printf("\n%d %d\n",bs[0][i][0],bs[0][i][1]);
    
  //printf("\n nn= %d , tdim = %d, nb =%d, nab =%d, nbs[0] = %d\n",npt,tdim,nb,nab,nbs[0]);
  
  //for (i = 0; i < npt; i++)
    //printf("\n %16.10e %16.10e \n",pt[i][0],pt[i][1]);

  int nt = -1;

  Point2D lo, hi;

  // determine extent of domain
  lo = Point2D(1.0e20,1.0e20);
  hi = Point2D(-1.0e20,-1.0e20);
  for (j=0; j < 2; j++)
  {
    for (n=0; n < npt; n++)
    {
      lo[j] = MIN(lo[j],pt[n][j]);
      hi[j] = MAX(hi[j],pt[n][j]);
    }
  }

  double ds;
  double tol = MIN(hi[0]-lo[0],hi[1]-lo[1])*1.0e-15;
  double dsmn = 1.0e20;
  for (b=0; b < nb; b++)
  {
    for (i=0; i < nbs[b]; i++)
    {
      n0 = bs[b][i][0];
      n1 = bs[b][i][1];
      ds = dist(pt[n0],pt[n1]);
      dsmn = MIN(ds,dsmn);
    }
  }
  if (dsmn < tol)
  {
    fprintf(stderr,"\nTRIMESH: zero length edge spacing encountered!");
    fprintf(stderr,"\n         tolerance value = %g",tol);
    fprintf(stderr,"\n         minimum edge spacing = %g\n",dsmn);
    for (b=0; b < nb; b++)
    {
      fprintf(stderr,"\n boundary %d:",b);
      for (i=0; i < nbs[b]; i++)
      {
        n0 = bs[b][i][0];
        n1 = bs[b][i][1];
        fprintf(stderr,"\n edge %d:",i);
        fprintf(stderr,"\n point %d",n0);
        pt[n0].print(stderr);
        fprintf(stderr,"\n point %d",n1);
        pt[n1].print(stderr);
      }
    }
    exit(0);
  }

  //dsmn = 1.0e20;
  //for (i=0; i < npt-1; i++)
  //{
  //  for (j=i+1; j < npt; j++)
  //  {
  //    ds = dist(pt[i],pt[j]);
  //    dsmn = MIN(ds,dsmn);
  //  }
  //}
  //if (dsmn < tol)
  //{
  //  fprintf(stderr,"\nTRIMESH: duplicate points encountered!");
  //  fprintf(stderr,"\n         tolerance value = %g",tol);
  //  fprintf(stderr,"\n         minimum distance = %g\n",dsmn);
  //  exit(0);
  //}

  // allocate space for local points
  Point2D *p;
  p = new Point2D[npt+4];

  // create bounding box
  double dx = (hi[0]-lo[0])*0.25;;
  double dy = (hi[1]-lo[1])*0.25;;
  lo -= Point2D(dx,dy);
  hi += Point2D(dx,dy);
  nn = 0;
  p[nn++] = Point2D(lo[0],lo[1]);
  p[nn++] = Point2D(hi[0],lo[1]);
  p[nn++] = Point2D(hi[0],hi[1]);
  p[nn++] = Point2D(lo[0],hi[1]);

  // create initial triangles
  nt=0;
  tri[nt][0] = 0;
  tri[nt][1] = 1;
  tri[nt][2] = 2;
  nt++;
  tri[nt][0] = 2;
  tri[nt][1] = 3;
  tri[nt][2] = 0;
  nt++;

  // initialize node to element hash table
  Lhash = new Linked_List*[npt+4];

  for (n=0; n < npt+4; n++)
    Lhash[n] = new Linked_List();

  mknbr = new Linked_List();

  for (t=0; t < nt; t++)
  {
    mknbr->Insert(t);
    for (i=0; i < 3; i++)
    {
      n = tri[t][i];
      Lhash[n]->Insert(t);
    }
  }

  int (*nbr)[3];
  nbr = new int[tdim][3];
  for (t=0; t < tdim; t++)
    nbr[t][0]=nbr[t][1]=nbr[t][2]=-1;

  make_nbrs(mknbr,nn,nt,tri,nbr,Lhash);

  Linked_List *del;

  int flag = 1;
  seed = 0;

  // insert points
  for (n=0; n < npt && flag; n++)
  {

    // find triangle containing point
    if (seed >= nt) seed = 0;
    seed = search(seed,pt[n],p,tri,nbr);

    assert(seed >= 0);

    m = seed;

    // create list of deleted triangles
    del = new Linked_List();

    //if (n == 12612)
    //{
    //  n0=tri[m][0];
    //  n1=tri[m][1];
    //  n2=tri[m][2];
    // FILE *fp;
    //  if ((fp = fopen("debug.dat","w")) == 0)
    //  {
    //    exit(0);
    //  }
    //  fprintf(fp,"%19.10e %19.10e 0.0\n",  p[n0][0],p[n0][1]);
    //  fprintf(fp,"%19.10e %19.10e 0.0\n",  p[n1][0],p[n1][1]);
    //  fprintf(fp,"%19.10e %19.10e 0.0\n",  p[n2][0],p[n2][1]);
    //  fprintf(fp,"%19.10e %19.10e 0.0\n\n",p[n0][0],p[n0][1]);
    //  fprintf(fp,"%19.10e %19.10e 0.0\n",  p[n0][0],p[n0][1]);
    //  fprintf(fp,"%19.10e %19.10e 0.0\n\n",pt[n][0],pt[n][1]);
    //  fprintf(fp,"%19.10e %19.10e 0.0\n",  p[n1][0],p[n1][1]);
    //  fprintf(fp,"%19.10e %19.10e 0.0\n\n",pt[n][0],pt[n][1]);
    //  fprintf(fp,"%19.10e %19.10e 0.0\n",  p[n2][0],p[n2][1]);
    //  fprintf(fp,"%19.10e %19.10e 0.0\n\n",pt[n][0],pt[n][1]);
    //  fclose(fp);
    //}

    // perform recursive neighbor search for other deleted triangles
    nbr_search(pt[n],p,tri,nbr,del,m);

    m = del->Length();

    if (m > 0)
    {
      // create list of saved edges
      int (*edge)[2];
      int edim = 3*m;
      edge = new int[edim][2];
      int ne = 0;

      hd = del->head;
      while (hd)
      {
        t = hd->data;
        for (i=0; i < 3; i++)
        {
          m = nbr[t][i];
          if (m < 0 || !del->In_list(m))
          {
            switch (i)
            {
              case 0: n0=tri[t][0]; n1=tri[t][1]; break;
              case 1: n0=tri[t][1]; n1=tri[t][2]; break;
              case 2: n0=tri[t][2]; n1=tri[t][0]; break;
            }
            assert(ne < edim);
            edge[ne][0] = n0;
            edge[ne][1] = n1;
            ne++;
          }
        }
        hd = hd->next;
      }

      // compress triangle list
      compress_list(del, Lhash, mknbr, tri, nbr, nt);

      // create new triangles
      p[nn] = pt[n];
      for (i=0; i < ne && flag; i++)
      {
        if (nt < tdim)
        {
          tri[nt][0] = nn;
          tri[nt][1] = edge[i][0];
          tri[nt][2] = edge[i][1];
          nbr[nt][0]=nbr[nt][1]=nbr[nt][2]=-1;

          // add new elements to node hash table
          for (j=0; j < 3; j++)
            Lhash[tri[nt][j]]->Insert(nt);

          mknbr->Insert(nt);
          nt++;
        } else
          flag = 0;
      }
      nn++;

      delete[] edge;
    } else
    {
      fprintf(stderr,"\nUh-Oh! Cannot find triangles to delete, likely due to winding issues!\n");
      fflush(stderr);
    }

    delete del;

    make_nbrs(mknbr,nn,nt,tri,nbr,Lhash);

  }

  //
  // write GNUPLOT file
  //
  //char filename[32];
  //FILE *fp;
  //filename[0]='\0';
  //strcat(filename,"trimesh1.dat");
  ////printf("\nFilename = <%s>",filename);
  //// Open file for write
  //if ((fp = fopen(filename,"w")) == 0)
  //{
  //  printf("\nError opening file <%s>.",filename);
  //  exit(0);
  //}
  //for (i=0; i < nt; i++)
  //{
  //  int n0 = tri[i][0];
  //  int n1 = tri[i][1];
  //  int n2 = tri[i][2];
  //  fprintf(fp,"%19.10e %19.10e 0.0\n",  p[n0][0],p[n0][1]);
  //  fprintf(fp,"%19.10e %19.10e 0.0\n",  p[n1][0],p[n1][1]);
  //  fprintf(fp,"%19.10e %19.10e 0.0\n",  p[n2][0],p[n2][1]);
  //  fprintf(fp,"%19.10e %19.10e 0.0\n\n",p[n0][0],p[n0][1]);
  //}
  //fclose(fp);

  //
  // ensure boundary edges & auxiliary edges are included
  //
  if (flag && (nb+nab) > 0)
  {
    Vector2D nrm, vec1;

    for (b=0; b < nb+nab && flag; b++)
    {
      for (int ib=0; ib < nbs[b] && flag; ib++)
      {
        n0 = bs[b][ib][0] + 4; // increment to account for superstructure nodes
        n1 = bs[b][ib][1] + 4; // increment to account for superstructure nodes
        Linked_Node *hd0, *hd1;
        hd0 = Lhash[n0]->head;
        m = 0;
        while (hd0 && !m)
        {
          t = hd0->data;
          hd1 = Lhash[n1]->head;
          while (hd1 && !m)
          {
            if (hd1->data == t)
              m = 1;
            hd1 = hd1->next;
          }
          hd0 = hd0->next;
        }
        if (!m) // edge not in mesh
        {
          Vector2D nrm1, nrm2, vs;
          vs = Vector2D(p[n0],p[n1]);
          vs.normalize();

          // walk through elements swapping to get to node n1
          n2 = -1;
          while (n2 != n1 && flag)
          {
            n2 = -1;
            hd0 = Lhash[n0]->head;
            while (hd0 && n2 < 0)
            {
              m = hd0->data;

              for (i=0; i < 3; i++)
              {
                switch (i)
                {
                  case 0:
                    m0=tri[m][0];
                    m1=tri[m][1];
                    m2=tri[m][2];
                    t=nbr[m][1];
                    break;
                  case 1:
                    m0=tri[m][1];
                    m1=tri[m][2];
                    m2=tri[m][0];
                    t=nbr[m][2];
                    break;
                  case 2:
                    m0=tri[m][2];
                    m1=tri[m][0];
                    m2=tri[m][1];
                    t=nbr[m][0];
                    break;
                }
                if (m0 == n0)
                  break;
              }

              // check for triangle containing vector VS
              nrm1 = Vector2D(p[m0][1]-p[m1][1],-(p[m0][0]-p[m1][0]));
              nrm1.normalize();
              nrm2 = Vector2D(p[m2][1]-p[m0][1],-(p[m2][0]-p[m0][0]));
              nrm2.normalize();
              if (t >= 0 && (vs * nrm1) > 1.0e-15 &&
                            (vs * nrm2) > 1.0e-15)
              {
                n2 = -1;
                for (i=0; i < 3 && n2 < 0; i++)
                {
                  switch (i)
                  {
                    case 0:
                      if (tri[t][0] == m2 && tri[t][1] == m1)
                        n2 = tri[t][2]; break;
                    case 1:
                      if (tri[t][1] == m2 && tri[t][2] == m1)
                        n2 = tri[t][0]; break;
                    case 2:
                      if (tri[t][2] == m2 && tri[t][0] == m1)
                        n2 = tri[t][1]; break;
                  }
                }
                assert(n2 >= 0);
            
                // reset neighbors and neighbor's neighbors
                for (i=0; i < 3; i++)
                {
                  n = nbr[m][i];
                  if (n >= 0)
                    mknbr->Insert(n);
                  n = nbr[t][i];
                  if (n >= 0)
                    mknbr->Insert(n);
                }
                mknbr->Insert(m);
                mknbr->Insert(t);
                    
                Lhash[m0]->Remove(m);
                Lhash[m1]->Remove(m);
                Lhash[m2]->Remove(m);
                Lhash[m1]->Remove(t);
                Lhash[m2]->Remove(t);
                Lhash[n2]->Remove(t);
                tri[m][0] = m0;
                tri[m][1] = m1;
                tri[m][2] = n2;
                Lhash[m0]->Insert(m);
                Lhash[m1]->Insert(m);
                Lhash[n2]->Insert(m);
                tri[t][0] = m0;
                tri[t][1] = n2;
                tri[t][2] = m2;
                Lhash[m0]->Insert(t);
                Lhash[n2]->Insert(t);
                Lhash[m2]->Insert(t);
                make_nbrs(mknbr,nn,nt,tri,nbr,Lhash);
              } else
                hd0 = hd0->next;
            }
            if (n2 < 0)
            {
              fprintf(stderr,"\nTrimesh: Edge not recovered!\n");
              fflush(stderr);
              flag = 0;
              //exit(0);
            }
          }
        }
      }
    }
  }

  //filename[0]='\0';
  //strcat(filename,"trimesh2.dat");
  //printf("\nFilename = <%s>",filename);
  // Open file for write
  //FILE *fp = NULL;
  //if ((fp = fopen("trimesh_swap.dat","w")) == 0)
  //{
    //fprintf(stderr,"\nError opening file <trimesh_swap.dat>.");
    //fflush(stderr);
    //exit(0);
  //}
  //for (i=0; i < nt; i++)
  //{
    //int n0 = tri[i][0];
    //int n1 = tri[i][1];
    //int n2 = tri[i][2];
    //fprintf(fp,"%19.10e %19.10e 0.0\n",  p[n0][0],p[n0][1]);
    //fprintf(fp,"%19.10e %19.10e 0.0\n",  p[n1][0],p[n1][1]);
    //fprintf(fp,"%19.10e %19.10e 0.0\n",  p[n2][0],p[n2][1]);
    //fprintf(fp,"%19.10e %19.10e 0.0\n\n",p[n0][0],p[n0][1]);
  //}
  //fclose(fp);
  
  //
  // delete exterior cells
  //
  if (flag)
  {
    // first delete any triangle attached to super-nodes (0-3)
    int *map = new int[nt];

    for (t=0; t < nt; t++)
      map[t] = 1;

    for (t=0; t < nt; t++)
    {
      n0 = tri[t][0];
      n1 = tri[t][1];
      n2 = tri[t][2];
      if (n0 < 4 || n1 < 4 || n2 < 4)
        map[t] = -1;
    }

    // compress triangle list
    for (m=t=0; t < nt; t++)
      if (map[t] > 0)
        map[t] = m++;
      else
        map[t] = -1;

    for (t=0; t < nt; t++)
    {
      if (map[t] >= 0)
      {
        tri[map[t]][0] = tri[t][0];
        tri[map[t]][1] = tri[t][1];
        tri[map[t]][2] = tri[t][2];
      }
    }
    nt = m;

    //filename[0]='\0';
    //strcat(filename,"trimesh3.dat");
    ////printf("\nFilename = <%s>",filename);
    //// Open file for write
    //if ((fp = fopen(filename,"w")) == 0)
    //{
    //  printf("\nError opening file <%s>.",filename);
    //  exit(0);
    //}
    //for (i=0; i < nt; i++)
    //{
    //  int n0 = tri[i][0];
    //  int n1 = tri[i][1];
    //  int n2 = tri[i][2];
    //  fprintf(fp,"%19.10e %19.10e 0.0\n",  p[n0][0],p[n0][1]);
    //  fprintf(fp,"%19.10e %19.10e 0.0\n",  p[n1][0],p[n1][1]);
    //  fprintf(fp,"%19.10e %19.10e 0.0\n",  p[n2][0],p[n2][1]);
    //  fprintf(fp,"%19.10e %19.10e 0.0\n\n",p[n0][0],p[n0][1]);
    //}
    //fclose(fp);

    if (nb > 0)
    {
      // recreate the node to triangle hash table
      for (n=0; n < npt+4; n++)
      {
        delete Lhash[n];
        Lhash[n] = new Linked_List();
      }
      for (t=0; t < nt; t++)
      {
        mknbr->Insert(t);
        nbr[t][0]=nbr[t][1]=nbr[t][2]=-1;
        for (i=0; i < 3; i++)
        {
          n = tri[t][i];
          Lhash[n]->Insert(t);
        }
      }

      make_nbrs(mknbr,nn,nt,tri,nbr,Lhash);

      // mark using painting algorithm
      for (t=0; t < nt; t++)
        map[t] = 0;

      for (b=0; b < nb; b++)
      {
        for (i=0; i < nbs[b]; i++)
        {
          m0 = bs[b][i][0] + 4; // increment to account for superstructure
          m1 = bs[b][i][1] + 4; // increment to account for superstructure

          // identify triangles containing both nodes
          Linked_Node *hd0, *hd1;
          hd0 = Lhash[m0]->head;
          while (hd0)
          {
            t = hd0->data;
            n0 = tri[t][0];
            n1 = tri[t][1];
            n2 = tri[t][2];
            hd1 = Lhash[m1]->head;
            //while (hd1 && map[t] == 0)
            while (hd1)
            {
              if (hd1->data == t)
              {
                // use node ordering to determine left/right(in/out) status
                if (((n0 == m0 && n1 == m1) ||
                    (n1 == m0 && n2 == m1) ||
                    (n2 == m0 && n0 == m1)) && map[t] == 0)
                  map[t] = 1;
                if ((n0 == m1 && n1 == m0) ||
                    (n1 == m1 && n2 == m0) ||
                    (n2 == m1 && n0 == m0))
                  map[t] = -1;
              }
              hd1 = hd1->next;
            }
            hd0 = hd0->next;
          }
        }
      }

      // perform flood-fill to notify neighbors
      do
      {
        for (t=0; t < nt; t++)
        {
          if (map[t] == 0)
            continue;
          for (i=0; i < 3; i++)
          {
            n = nbr[t][i];
            if (n < 0)
              continue;
            if (map[n] == 0)
              map[n] = map[t];
          }
        }
        for (m=t=0; t < nt; t++)
          if (map[t] == 0)
            m++;
      } while (m > 0);
        
      // compress triangle list
      for (m=t=0; t < nt; t++)
        if (map[t] > 0)
          map[t] = m++;

      for (t=0; t < nt; t++)
      {
        if (map[t] >= 0)
        {
          tri[map[t]][0] = tri[t][0];
          tri[map[t]][1] = tri[t][1];
          tri[map[t]][2] = tri[t][2];
        }
      }
      nt = m;
    }

    delete[] map;
  }

  //
  // write GNUPLOT file
  //
  //char filename[32];
  //FILE *fp;
  //filename[0]='\0';
  //strcat(filename,"trimesh4.dat");
  ////printf("\nFilename = <%s>",filename);
  //// Open file for write
  //if ((fp = fopen(filename,"w")) == 0)
  //{
  //  printf("\nError opening file <%s>.",filename);
  //  exit(0);
  //}
  //for (i=0; i < nt; i++)
  //{
  //  int n0 = tri[i][0];
  //  int n1 = tri[i][1];
  //  int n2 = tri[i][2];
  //  fprintf(fp,"%19.10e %19.10e 0.0\n",  p[n0][0],p[n0][1]);
  //  fprintf(fp,"%19.10e %19.10e 0.0\n",  p[n1][0],p[n1][1]);
  //  fprintf(fp,"%19.10e %19.10e 0.0\n",  p[n2][0],p[n2][1]);
  //  fprintf(fp,"%19.10e %19.10e 0.0\n\n",p[n0][0],p[n0][1]);
  //}
  //fclose(fp);

  delete[] nbr;
  delete mknbr;

  // decrement all node indices by 4 to account for superstructure
  for (t=0; t < nt; t++)
    for (i=0; i < 3; i++)
      tri[t][i] -= 4;

  if (flag)
    flag = nt;
  else
    flag = -1;

  for (n=0; n < npt+4; n++)
    delete Lhash[n];
  delete[] Lhash;

  delete[] p;

  //
  // write GNUPLOT file
  //
  //char filename[32];
  //filename[0]='\0';
  //strcat(filename,"trimesh.dat");
  ////printf("\nFilename = <%s>",filename);
  //// Open file for write
  // FILE *fp;
  //if ((fp = fopen(filename,"w")) == 0)
  //{
  //  printf("\nError opening file <%s>.",filename);
  //  exit(0);
  //}
  //for (i=0; i < nt; i++)
  //{
  //  int n0 = tri[i][0];
  //  int n1 = tri[i][1];
  //  int n2 = tri[i][2];
  //  fprintf(fp,"%19.10e %19.10e 0.0\n",  pt[n0][0],pt[n0][1]);
  //  fprintf(fp,"%19.10e %19.10e 0.0\n",  pt[n1][0],pt[n1][1]);
  //  fprintf(fp,"%19.10e %19.10e 0.0\n",  pt[n2][0],pt[n2][1]);
  //  fprintf(fp,"%19.10e %19.10e 0.0\n\n",pt[n0][0],pt[n0][1]);
  //}
  //fclose(fp);

  return(flag);
}

void compress_list(Linked_List *del, Linked_List **Lhash, Linked_List *mknbr,
                   int tri[][3], int nbr[][3], int &nt)
{
  int i, j, n, t;

  // compress triangle list
  Linked_Node *hd;

  hd = del->head;
  while (hd)
  {
    t = hd->data;
    mknbr->Remove(t);

    // delete current element from neighbor's connectivity
    for (i=0; i < 3; i++)
    {
      Lhash[tri[t][i]]->Remove(t);
      n = nbr[t][i];
      if (n >= 0)
      {
        if (!del->In_list(n) && !mknbr->In_list(n))
          mknbr->Insert(n);
        for (j=0; j < 3; j++)
          if (nbr[n][j] == t)
            nbr[n][j] = -1;
      }
    }
    hd = hd->next;
  }

  while ((hd = del->head))
  {
    t = hd->data;

    while (del->In_list(nt-1))
    {
      mknbr->Remove(nt-1);

      // delete current last element
      del->Remove(nt-1);
      nt--;

      // reset connectivities
      tri[nt][0]=tri[nt][1]=tri[nt][2]=-1;
      nbr[nt][0]=nbr[nt][1]=nbr[nt][2]=-1;
    }

    if (t < nt-1)
    {
      
      // change neighbor's connectivity
      if (mknbr->In_list(nt-1))
        mknbr->Replace(nt-1,t);
      for (i=0; i < 3; i++)
      {
        Lhash[tri[nt-1][i]]->Remove(nt-1);
        n = nbr[nt-1][i];
        if (n >= 0)
        {
          if (!del->In_list(n) && !mknbr->In_list(n))
            mknbr->Insert(n);
          for (j=0; j < 3; j++)
            if (nbr[n][j] == nt-1)
              nbr[n][j] = t;
        }
      }
      // transfer nodes and neighbors
      tri[t][0]=tri[nt-1][0];
      tri[t][1]=tri[nt-1][1];
      tri[t][2]=tri[nt-1][2];
      nbr[t][0]=nbr[nt-1][0];
      nbr[t][1]=nbr[nt-1][1];
      nbr[t][2]=nbr[nt-1][2];
      if (!mknbr->In_list(t))
        mknbr->Insert(t);
      for (i=0; i < 3; i++)
      {
        Lhash[tri[t][i]]->Insert(t);
      }
      // delete current element from delete list
      del->Remove(t);
      nt--;
      // reset connectivities
      tri[nt][0]=tri[nt][1]=tri[nt][2]=-1;
      nbr[nt][0]=nbr[nt][1]=nbr[nt][2]=-1;
    }

  }
}

void nbr_search(Point2D pt, Point2D p[], int tri[][3], int nbr[][3], Linked_List *del, int m)
{
  if (del->In_list(m))
    return;

  int n0 = tri[m][0];
  int n1 = tri[m][1];
  int n2 = tri[m][2];
  if (circle_test(p[n0],p[n1],p[n2],pt))
  {
    del->Insert(m);
    if (nbr[m][0] >= 0)
      nbr_search(pt,p,tri,nbr,del,nbr[m][0]);
    if (nbr[m][1] >= 0)
      nbr_search(pt,p,tri,nbr,del,nbr[m][1]);
    if (nbr[m][2] >= 0)
      nbr_search(pt,p,tri,nbr,del,nbr[m][2]);
  }

  return;
}

void make_nbrs(Linked_List *mknbr, int nn, int nt, int tri[][3], int nbr[][3], Linked_List **Lhash)
{
  Linked_Node *hd;
  int c, m, n, s;
  int n0, n1;

  while ((hd = mknbr->head))
  {
    c = hd->data;

    for (s=0; s < 3; s++)
    {
      switch (s)
      {
        case 0: n0 = tri[c][0]; n1 = tri[c][1]; break;
        case 1: n0 = tri[c][1]; n1 = tri[c][2]; break;
        case 2: n0 = tri[c][2]; n1 = tri[c][0]; break;
      }
      m = -1;
      Linked_Node *hd0, *hd1;
      hd0 = Lhash[n0]->head;
      while (hd0 && m < 0)
      {
        n = hd0->data;
        if (n != c)
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
      nbr[c][s] = m;
    }
    mknbr->Remove(c);
  }

}

double det_3X3(double m[3][3])
{
  double det;

  det = m[0][0]*m[1][1]*m[2][2]+
        m[0][1]*m[1][2]*m[2][0]+
        m[0][2]*m[1][0]*m[2][1]-
        m[0][2]*m[1][1]*m[2][0]-
        m[0][1]*m[1][0]*m[2][2]-
        m[0][0]*m[1][2]*m[2][1];

  return(det);
}

int circle_test(Point2D t1, Point2D t2, Point2D t3, Point2D t)
{
  int flag = 0;
  Point2D p1, p2, p3, p;

  // scale incoming triangle to unit max length
  Point2D lo, hi;
  lo[0] = MIN(t1[0],MIN(t2[0],t3[0]));
  lo[1] = MIN(t1[1],MIN(t2[1],t3[1]));
  hi[0] = MAX(t1[0],MAX(t2[0],t3[0]));
  hi[1] = MAX(t1[1],MAX(t2[1],t3[1]));
  double dm = MAX(hi[0]-lo[0],hi[1]-lo[1]);
  dm = MAX(1.0e-15,dm);
  p1[0] = (t1[0]-lo[0])/dm;
  p1[1] = (t1[1]-lo[1])/dm;
  p2[0] = (t2[0]-lo[0])/dm;
  p2[1] = (t2[1]-lo[1])/dm;
  p3[0] = (t3[0]-lo[0])/dm;
  p3[1] = (t3[1]-lo[1])/dm;
  p[0] = (t[0]-lo[0])/dm;
  p[1] = (t[1]-lo[1])/dm;

//
//  compute circumscribing radius
//
  Point2D cc;

  double m[3][3];
  double a, bx, by, c;

  m[0][0] = p1[0];
  m[1][0] = p2[0];
  m[2][0] = p3[0];
  m[0][1] = p1[1];
  m[1][1] = p2[1];
  m[2][1] = p3[1];
  m[0][2] = 1.0;
  m[1][2] = 1.0;
  m[2][2] = 1.0;
  a = det_3X3(m);
  a = SIGN(MAX(1.0e-15,fabs(a)),a);

  m[0][0] = p1[0]*p1[0]+p1[1]*p1[1];
  m[1][0] = p2[0]*p2[0]+p2[1]*p2[1];
  m[2][0] = p3[0]*p3[0]+p3[1]*p3[1];
  bx = det_3X3(m);

  m[0][1] = p1[0];
  m[1][1] = p2[0];
  m[2][1] = p3[0];
  by = det_3X3(m);

  m[0][2] = p1[1];
  m[1][2] = p2[1];
  m[2][2] = p3[1];
  c = det_3X3(m);

  cc = Point2D(0.5*bx/a,-0.5*by/a);

  double rc, rp;
  Vector2D vc, vp;
  vp = Vector2D(cc,p);
  rp = vp.magnitude();
  //vc = Vector2D(cc,p1);
  //rc = vc.magnitude();
  //vc = Vector2D(cc,p2);
  //rc = vc.magnitude();
  //vc = Vector2D(cc,p3);
  //rc = vc.magnitude();

  rc = 0.5*sqrt(bx*bx+by*by+4.0*a*c)/fabs(a);

  //if (rp <= 1.00001*rc)
  if (rp <= rc)
    flag = 1;

  return (flag);
}

double dist( Point2D p1, Point2D p2)
{
  double dx = p2[0]-p1[0];
  double dy = p2[1]-p1[1];
  double mag = sqrt(dx*dx+dy*dy);
  return mag;
}

int search(int t, Point2D pt, Point2D p[], int tri[][3], int nbr[][3])
{
  int n, n0, n1, n2;
  Vector2D v1, v2, v3, nrm1, nrm2, nrm3;
  //double tol = 1.0e-20;

  n0 = tri[t][0];
  n1 = tri[t][1];
  n2 = tri[t][2];
  nrm1 = Vector2D((p[n1][1]-p[n0][1]),-(p[n1][0]-p[n0][0]));
  nrm2 = Vector2D((p[n2][1]-p[n1][1]),-(p[n2][0]-p[n1][0]));
  nrm3 = Vector2D((p[n0][1]-p[n2][1]),-(p[n0][0]-p[n2][0]));
  nrm1.normalize();
  nrm2.normalize();
  nrm3.normalize();
  v1 = Vector2D((p[n0]+p[n1])*0.5,pt);
  v2 = Vector2D((p[n1]+p[n2])*0.5,pt);
  v3 = Vector2D((p[n2]+p[n0])*0.5,pt);
  v1.normalize();
  v2.normalize();
  v3.normalize();
  double dot1 = v1*nrm1;
  double dot2 = v2*nrm2;
  double dot3 = v3*nrm3;
  if (dot1 <= 0.0 && dot2 <= 0.0 && dot3 <= 0.0)
    n = t;
  else
  {
    if (nbr[t][0] >= 0 && dot1 >= MAX(dot2,dot3))
      n = search(nbr[t][0],pt,p,tri,nbr);
    else if (nbr[t][1] >= 0 && dot2 >= MAX(dot1,dot3))
      n = search(nbr[t][1],pt,p,tri,nbr);
    else if (nbr[t][2] >= 0 && dot3 >= MAX(dot1,dot2))
      n = search(nbr[t][2],pt,p,tri,nbr);
    else
      n = -1;
  }

  return(n);
}

double max_angle( Point2D p1, Point2D p2, Point2D p3)
{
  double dot, mag;
  Vector2D v1, v2;

  v1 = Vector2D(p1,p2);
  v2 = Vector2D(p1,p3);
  v1.normalize();
  v2.normalize();
  dot = v1 * v2;
  mag = acos(dot);

  v1 = Vector2D(p2,p3);
  v2 = Vector2D(p2,p1);
  v1.normalize();
  v2.normalize();
  dot = v1 * v2;
  mag = MAX(mag,acos(dot));

  v1 = Vector2D(p3,p1);
  v2 = Vector2D(p3,p2);
  v1.normalize();
  v2.normalize();
  dot = v1 * v2;
  mag = MAX(mag,acos(dot));

  return mag;
}
#ifdef SPACE
#undef SPACE
#endif
