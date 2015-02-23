#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include "geometry.h"
#include "FASTAR.h"
#include "smooth.h"
#include "Util.h"
#include "trimesh.h"
#include "Spacing_Field.h"
#include "PList.h"
#include "Octree_Storage.h"
//#include "tetgen.h"

//global variables
extern int my_rank; /*rank of process*/
extern int num_procs; /*number of processes*/
extern int ugeom; //spacing file only (1) or both spacing and geom used
extern int adapt; //using adaptation or not
extern int restart; //restarting or not
extern int ftype; //0 for UCD, 1 for STL
extern int ugeom; //spacing file only (1) or both spacing and geom used
extern int evencut; //if we want nice box around geometry
extern int bb; //num bds to be cut
extern int gentet; //use TetGen
extern int viscous; //insert layersextern 
extern double factor; //geom prog factor
extern double vspace; //viscous spacing
extern int nvbd; //number viscous bd
extern int highquality; //allows desirable constraints to be enforced
extern int numregions; //gives number of regions to keep from tetgen
extern int overset; //overset or not
extern FILE *in_f, *jou_f, *out_f, *debug_f; /*input output journal files global*/

void fastar_lib(geometry *geom, double smn, double smx, double ar, double ctol, double mtol, double nspace, int *boxcut, int *vbd, int *regions)
{
  int bdim = 132;
  char *sname = new char[bdim];
  int i, j, k, l, m, n;
  FILE *ucd_f = NULL;
  FILE *tet_f = NULL;
  //declared here for use w or w/out restart
  POLYMESH *finalmesh = new POLYMESH();
  int n0, n1, n2, snn = 0, gnn = 0, nni = 0, flag = 0, ngt = 0, tempnt = 0;
  double tol, gtol;
  const int bdim2 = 400;
  char buff[bdim2];
  Vector v0, v1, norm, normv;
  int *map = NULL;
  Point pv, p0, p1, cent, fcent;
  double olo[3], ohi[3];
  time_t tm;
  char *t_char;
  int nonviscfacets = 0;
  //created and deleted regardless
  POLYMESH *ovst = new POLYMESH();

  if (restart == 0)
  {
  if (my_rank == 0)
  {
    fprintf(out_f,"\nBeginning mesh generation.\n");
    fflush(out_f);
  }
  
  //need to declare stmesh
  split_tree *otmesh = new split_tree();

  time(&tm);
  t_char = ctime(&tm);
  fprintf(out_f,"\nBegan tree creation at %s",t_char);
    
  //create tree
  otmesh->create_tree(geom, smn, smx, ar, mtol, nspace);

  time(&tm);
  t_char = ctime(&tm);
  fprintf(out_f,"\nFinished tree creation and started cutting at %s",t_char);

  //cut mesh
  if(!overset)
    otmesh->cut_voxels(geom, ctol, smn, smx, boxcut, vbd);

  time(&tm);
  t_char = ctime(&tm);
  fprintf(out_f,"\nFinished cutting at %s",t_char);

  //if overset, mimic cutting
  for (i = 0; i < otmesh->nvoxels && overset; i++)
    if (otmesh->voxels[i].cut != 2)
      otmesh->voxels[i].cut = 1;

  if (my_rank == 0)
  {
    fprintf(out_f,"\nFinished mesh generation, outputting volume and surface meshes....\n");
    fflush(out_f);
  }
  
  //finally, output mesh
  
  //init values
  finalmesh->nn = 0;
  finalmesh->nelems = 0;
  finalmesh->nb = 0;

  ovst->nn = 0;
  ovst->nelems = 0;
  ovst->nb = 0;

#ifdef _NONUNIQUE
  Point *pt = new Point[20];
  double begx, begy, begz, endx, endy, endz, halfx, halfy, halfz;
  Point p0, p1;
  int nde[4];

  //NON-UNIQUE, but only finest level
  for (k = 0; k < otmesh->nvoxels; k++)
    {
    if (otmesh->voxels[k].children.max > 0 || otmesh->voxels[k].cut != 1)
      continue;
    
    finalmesh->nelems++; //set up to init
      
    //set all nodes
    for (i = 0; i < 8; i++)
      {
      otmesh->voxels[k].nodes[i] = finalmesh->nn;
      finalmesh->nn++;
      }
    }
    
  fprintf(out_f,"\nNumber of elements = %d, number of voxels = %d, number of nodes = %d\n",finalmesh->nelems,otmesh->nvoxels,finalmesh->nn);
  fflush(out_f);
    
  //now, set up element array
  finalmesh->element = new POLY_ELEMENT[finalmesh->nelems];
  //now, set up node array
  finalmesh->node = new NODE[finalmesh->nn];
    
  for (k = 0; k < otmesh->nvoxels; k++)
    {
    if (otmesh->voxels[k].children.max > 0 || otmesh->voxels[k].cut != 1)
      continue;
    
    p0 = otmesh->voxels[k].cornerpts[0];
    p1 = otmesh->voxels[k].cornerpts[1];
    
    //set up values for orig cornerpts
    begx = p0[0];
    begy = p0[1];
    begz = p0[2];
    endx = p1[0];
    endy = p1[1];
    endz = p1[2];
    halfx = (p0[0]+p1[0])/2.0;
    halfy = (p0[1]+p1[1])/2.0; 
    halfz = (p0[2]+p1[2])/2.0;
          
    //set up other pts of orig voxel
    pt[0] = p0;
    pt[1] = Point(endx,begy,begz); 
    pt[2] = Point(endx,endy,begz);
    pt[3] = Point(begx,endy,begz);
    pt[4] = Point(begx,begy,endz);
    pt[5] = Point(endx,begy,endz);
    pt[6] = p1; 
    pt[7] = Point(begx,endy,endz);
    pt[8] = Point(halfx,begy,begz);
    pt[9] = Point(endx,halfy,begz);
    pt[10] = Point(halfx,endy,begz);
    pt[11] = Point(begx,halfy,begz);
    pt[12] = Point(begx,begy,halfz);
    pt[13] = Point(endx,begy,halfz);
    pt[14] = Point(endx,endy,halfz);
    pt[15] = Point(begx,endy,halfz);
    pt[16] = Point(halfx,begy,endz);
    pt[17] = Point(endx,halfy,endz);
    pt[18] = Point(halfx,endy,endz);
    pt[19] = Point(begx,halfy,endz);
      
    //set all nodes
    for (i = 0; i < 8; i++)
      {
      finalmesh->node[otmesh->voxels[k].nodes[i]].vert = pt[i];
      }
    }
    
  delete [] pt;
    
  j = 0;  
  //now, set up elements
  for (i = 0; i < otmesh->nvoxels; i++)
    {
    if (otmesh->voxels[i].cut != 1 || otmesh->voxels[i].children.max > 0)
      continue;
      
    //now, set up proper lists
    finalmesh->element[j].nf = 6;
    finalmesh->element[j].f_n = new List*[6];
    for (k = 0; k < 6; k++)
      finalmesh->element[j].f_n[k] = new List();
      
    for (k = 0; k < 6; k++)
      {
      switch (k) 
        {
        case 0:
          nde[0] = 0;
          nde[1] = 3;
          nde[2] = 2;
          nde[3] = 1;
        break;
        case 1:
          nde[0] = 0;
          nde[1] = 1;
          nde[2] = 5;
          nde[3] = 4;
        break;
        case 2:
          nde[0] = 1;
          nde[1] = 2;
          nde[2] = 6;
          nde[3] = 5;
        break;
        case 3:
          nde[0] = 2;
          nde[1] = 3;
          nde[2] = 7;
          nde[3] = 6;
        break;
        case 4:
          nde[0] = 3;
          nde[1] = 0;
          nde[2] = 4;
          nde[3] = 7;
        break;
        case 5:
          nde[0] = 6;
          nde[1] = 7;
          nde[2] = 4;
          nde[3] = 5;
        break;
        default:
        break;
        }
        
      for (n = 0; n < 4; n++)  
        finalmesh->element[j].f_n[k]->Add_To_List(otmesh->voxels[i].nodes[nde[n]]);
      }
    j++;
    }    
  } //for restart == 0
#else
  Point *pts = new Point[8];
  int *nums = new int[8];
  int tnn = 0;

  time(&tm);
  t_char = ctime(&tm);
  fprintf(out_f,"\nStarted distinct node creation at %s",t_char);

  //need to count elements, generate nodes/physical nodes, set up elements
  for (k = 0; k < otmesh->nvoxels; k++)
    {
    if (otmesh->voxels[k].children.max == 0 && otmesh->voxels[k].cut == 1)
      finalmesh->nelems++; //set up to init
      
    //init all nodes
    for (i = 0; i < 27; i++)
      otmesh->voxels[k].nodes[i] = -1;
    }
  
  int tempne = 0;
  //now, set up element array...need to use my_mem, so we can add tets
  my_mem(finalmesh->element_dim, tempne, &(finalmesh->element), finalmesh->nelems, ELEMENT_CHUNK);
  for (k = tempne; k < finalmesh->nelems; k++)
    finalmesh->element[k].Initialize();
    
  //eventually will need to add in tet elements from cut region...cut = 2, uncut in = 1, uncut out = -1, undecided = 0
  fprintf(out_f,"\nNumber of elements = %d, number of voxels = %d\n",finalmesh->nelems,otmesh->nvoxels);
  fflush(out_f);
  
  //now, set up nodes only on finest level
  for (k = 0; k < otmesh->nvoxels; k++)
    {
    if (otmesh->voxels[k].children.max == 0 && otmesh->voxels[k].cut > 0)
      {
      //init nums
      for (i = 0; i < 8; i++)
        nums[i] = -1;
        
      //fprintf(out_f, "\nvoxel %d (%lf, %lf, %lf) (%lf, %lf, %lf)\n",k,otmesh->voxels[k].cornerpts[0][0],otmesh->voxels[k].cornerpts[0][1],otmesh->voxels[k].cornerpts[0][2],otmesh->voxels[k].cornerpts[1][0],otmesh->voxels[k].cornerpts[1][1],otmesh->voxels[k].cornerpts[1][2]);
      //fflush(out_f);
          
      //acutally set nodes
      otmesh->set_nodes(k, smn, smx, tnn, nums, pts);
      
      //now, need to create nodes
      if (tnn > finalmesh->node_dim)
        {
        //allocate mem
        my_mem(finalmesh->node_dim, finalmesh->nn, &(finalmesh->node), tnn, NODE_CHUNK);
        }
        
      for (i = 0; i < 8; i++)
        {
        if (nums[i] < 0)
          continue;
          
        //create physical node
        finalmesh->node[nums[i]].vert = pts[i];
        }
        
      //finally, set finalmesh->nn
      finalmesh->nn = tnn;
      } 
    }
    
  //since setting up nodes on cut and in (since cut sometimes finer level), need to ditch unused nodes from cut interface with cut/out
  int *map1 = new int[finalmesh->nn];
  
  //set as -1
  for (i = 0; i < finalmesh->nn; i++)
    map1[i] = -1;

  //set up map
  for (i = 0; i < otmesh->nvoxels; i++)
    if (otmesh->voxels[i].cut == 1)
      for (j = 0; j < 27; j++)
        if (otmesh->voxels[i].nodes[j] >= 0)
          map1[otmesh->voxels[i].nodes[j]] = 1;
   
  //set new counter
  int num = 0;
                 
  //now, reset nodes
  for (i = 0; i < finalmesh->nn; i++)
    if (map1[i] == 1)
      map1[i] = num++;
      
  //now, reset voxels nodes
  for (i = 0; i < otmesh->nvoxels; i++)
    if (otmesh->voxels[i].cut == 1)
      for (j = 0; j < 27; j++)
        if (otmesh->voxels[i].nodes[j] >= 0)
          otmesh->voxels[i].nodes[j] = map1[otmesh->voxels[i].nodes[j]];
          
  //finally, reset physical nodes, and realloc
  for (i = 0; i < finalmesh->nn; i++)
    {
    if (map1[i] >= 0 && map1[i] != i)
      {
      finalmesh->node[map1[i]].vert = finalmesh->node[i].vert;
      }
    }
  //fix mem  
  if (num < finalmesh->nn)
     my_mem(finalmesh->node_dim, finalmesh->nn, &(finalmesh->node), num, NODE_CHUNK);
  //reset nn
  finalmesh->nn = num;
  
  delete [] map1;
  
  fprintf(out_f,"\nNumber of nodes = %d\n",finalmesh->nn);
  fflush(out_f);
  
  delete [] pts;
  delete [] nums;

  time(&tm);
  t_char = ctime(&tm);
  fprintf(out_f,"\nFinished distinct node creation at %s",t_char);

  fprintf(out_f,"\nAllocating for geometry and stitching boundary objects.\n");
  fflush(out_f);
  
  //set number of boundaries...make last one internal bd, go back and delete later
  if (overset)
    {
    finalmesh->nb = 1;
    ovst->nb = nvbd+1;
    }
  else if (viscous && !overset)
    finalmesh->nb = geom->ngb+2;
  else
    finalmesh->nb = geom->ngb+1; 
  
  //need to add boundaries
  if (overset)
    {
    finalmesh->boundary = new BOUNDARY_MESH[1];
    ovst->boundary = new BOUNDARY_MESH[nvbd+1];
    }
  else if (viscous && !overset)
    finalmesh->boundary = new BOUNDARY_MESH[geom->ngb + 2];
  else
    finalmesh->boundary = new BOUNDARY_MESH[geom->ngb + 1];
    
  //first, set all names
  k = 0;
  for (i = 0; i < geom->ngb; i++)
    {
    if (overset)
      {
      flag = 0;
      for (j = 0; j < nvbd && !flag; j++)
        if (i+1 == vbd[j])
          flag++;
      }
    if (overset && flag)
      {
      ovst->boundary[k].name = new char[400];
      sprintf(ovst->boundary[k].name,"%s",geom->g_bname[i+1]);
      k++;
      }
    else
      {
      finalmesh->boundary[i].name = new char[400];
      sprintf(finalmesh->boundary[i].name,"%s",geom->g_bname[i+1]);
      }
    }

  if (!overset)
    { 
    finalmesh->boundary[geom->ngb].name = new char[400];
    sprintf(finalmesh->boundary[geom->ngb].name,"Stitching Boundary");
    }
  else
    { 
    ovst->boundary[0].name = new char[400];
    sprintf(ovst->boundary[0].name,"Outer Boundary");
    }

  if (viscous && !overset)
    {
    finalmesh->boundary[geom->ngb+1].name = new char[400];
    sprintf(finalmesh->boundary[geom->ngb+1].name,"Viscous Stitching Boundary");
    }
  if (viscous && overset)
    {
    ovst->boundary[nvbd].name = new char[400];
    sprintf(ovst->boundary[nvbd].name,"Viscous Stitching Boundary");
    }
     
  //now, set up lists!
  if (viscous && overset)
    {
    finalmesh->boundary[0].elist = new List();
    for (i = 0; i < nvbd+1; i++)
      {
      ovst->boundary[i].elist = new List();
      }
    }
  else if (viscous && !overset)
    {
    for (i = 0; i < geom->ngb+2; i++)
      {
      finalmesh->boundary[i].elist = new List();
      }
    }
  else
    {
    for (i = 0; i < geom->ngb+1; i++)
      {
      finalmesh->boundary[i].elist = new List();
      }
    }
      
  fprintf(out_f,"\nFinished allocating for geometry and stitching boundary objects.\n");
  fflush(out_f);

  int tnt = 0;
  List *tnode = new List();
  //be sure tempnt zeroed
  tempnt = 0;

  List nonviscbd;
  nonviscbd.construct();

  List keepfacet;
  keepfacet.construct();
 
  if (viscous)
   {
   time(&tm);
   t_char = ctime(&tm);
   fprintf(out_f,"\nStarted viscous extrusion creation at %s",t_char);
   //now, set up element array...need to use my_mem, so we can add tets
   //this will overalloc if not all bd viscous, but better more than not enough

   //now, set up node array
   if (!overset)
     my_mem(finalmesh->node_dim, finalmesh->nn, &(finalmesh->node), finalmesh->nn+(viscous*geom->n_gnodes), NODE_CHUNK);
   else
     my_mem(ovst->node_dim, ovst->nn, &(ovst->node), viscous*geom->n_gnodes, NODE_CHUNK);

   //now, create stack lists
   List **stack = new List*[geom->n_gnodes];
   for (i = 0; i < geom->n_gnodes; i++)
     stack[i] = new List();

   //we need a list of shared facets
   //since this just works for convex, we use 45 degrees as rotation
   //we need to be sure this is observed on all bds, not just viscous
   List shared;
   shared.construct();

   List bdnodes;
   bdnodes.construct();

   if (!overset)
   {
   for (i = 1; i <= geom->ngb+2; i++)
     {
     flag = 0; 
     for (j = 0; j < nvbd && !flag; j++)
       {
       if (i == vbd[j])
         flag++;
       }
     if (!flag)
       {
       nonviscbd.Add_To_List(i);
       }
     }
   }
   else
   {
   for (i = 1; i <= geom->ngb; i++)
     {
     flag = 0; 
     for (j = 0; j < nvbd && !flag; j++)
       {
       if (i == vbd[j])
         flag++;
       }
     if (!flag)
       {
       nonviscbd.Add_To_List(i);
       }
     }
   }

   nonviscbd.print(out_f);

   double nstore[3];

   int **tempv = new int*[geom->n_gnodes];
   for (i = 0; i < geom->n_gnodes; i++)
     {
     tempv[i] = new int[geom->ngb];
     for (j = 0; j < geom->ngb; j++)
       tempv[i][j] = -1;
     }

   for (i = 1; i <= geom->ngb; i++)
     for (j=geom->n_begin[i]; j < geom->n_end[i]; j++)
       for (k = 0; k < 3; k++)
         tempv[geom->g_facet[j].nodes[k]][i-1] = j;  //gives last facet visited, but good enough for norm avg

   int *omap = new int [geom->n_gnodes];

   //init
   for (i = 0; i < geom->n_gnodes && overset; i++)
     omap[i] = -1; 

   int flagv = 0;
   int new_node = 0;
   for (i = 0; i < geom->n_gnodes; i++)
     {
     flagv = 0;
     for (j = 0; j < geom->ngb; j++)
       {
       if (tempv[i][j] >= 0)
         flagv++;
       if (overset && tempv[i][j] >= 0 && omap[i] < 0 && !nonviscbd.Is_In_List(j+1))
         omap[i] = new_node++;
       }
     if (flagv > 1)
       shared.Add_To_List(i);
     }

   Vector normt;

   //reset ovst nn
   ovst->nn = new_node;

   fprintf(out_f,"\novst->nn = %d\n",ovst->nn);
   fflush(out_f);
 
   //if overset, create bd elem, remap bd nodes and add in
   for (i = 0; i < geom->n_gnodes && overset; i++)
     {
     if (omap[i] < 0)
       continue;
     ovst->node[omap[i]].vert = geom->g_vert[i];
     }

   if (overset && ((ovst->nelems+geom->n_gfacets) > ovst->element_dim))
    my_mem(ovst->element_dim, ovst->nelems, &(ovst->element), (ovst->nelems+geom->n_gfacets), ELEMENT_CHUNK);

   //since all elements are one face, set up nf and f_n now
   for (k = ovst->nelems; k < (ovst->nelems+geom->n_gfacets) && overset; k++)
      {
      ovst->element[k].Initialize();
      }

   k= 0;
   for (i = 1; i <= geom->ngb && overset; i++)
     {
     if (nonviscbd.Is_In_List(i))
       continue;
   
     for (j=geom->n_begin[i]; j < geom->n_end[i]; j++)
       {
       n0 = omap[geom->g_facet[j].nodes[0]];
       n1 = omap[geom->g_facet[j].nodes[1]];
       n2 = omap[geom->g_facet[j].nodes[2]];

       ovst->element[ovst->nelems].nf = 1;
       ovst->element[ovst->nelems].f_n = new List*[1];
       ovst->element[ovst->nelems].f_n[0] = new List();

       //take the resulting triangles and set up tri
       ovst->element[ovst->nelems].f_n[0]->Add_To_List(n0);
       ovst->element[ovst->nelems].f_n[0]->Add_To_List(n1);
       ovst->element[ovst->nelems].f_n[0]->Add_To_List(n2);

       ovst->boundary[k].elist->Add_To_List(ovst->nelems);
       ovst->nelems++;
       }
     k++;
     }

   int vpt1 = 0, vpt2 = 0, vpt3 = 0, bout = 0, tst = 0;

   for (i = 0; i < nvbd; i++)
     {
     for (j=geom->n_begin[vbd[i]]; j < geom->n_end[vbd[i]]; j++)
       {
       //reset bd check
       vpt1 = vpt2 = vpt3 = bout = 0;

       n0 = geom->g_facet[j].nodes[0];
       n1 = geom->g_facet[j].nodes[1];
       n2 = geom->g_facet[j].nodes[2];

       //if we border a nonvisc bd, we will make a special element
       vpt1 = shared.Is_In_List(n0);
       vpt2 = shared.Is_In_List(n1);
       vpt3 = shared.Is_In_List(n2);

       if (vpt1 || vpt2 || vpt3)
         {
         for (m = 0; m < 3 && !bout; m++)
           {
           tst = 0;
           switch(m)
             {
             case 0:
               n = n0;
               tst = vpt1;
             break;
             case 1:
               n = n1;
               tst = vpt2;
             break;
             case 2:
               n = n2;
               tst = vpt3;
             break;
             default:
             break;
             }
           if (tst)
             {
             for (k = 0; k < geom->ngb && !bout; k++)
               {
               if (tempv[n][k] >= 0 && nonviscbd.Is_In_List(k+1))
                 bout++;
               }
             }
           }
         }

       //if !bout, then shared amongst viscous
       if (bout)
         {
         //here, we add facet to list, so when we add geom pieces later, it is not skipped, added to right list, and pts added properly
         keepfacet.Add_To_List(j);
         //no need to create stack
         continue;  
         }

       v0 = Vector(geom->g_vert[n0],geom->g_vert[n1]);
       v1 = Vector(geom->g_vert[n0],geom->g_vert[n2]);
       norm = v0 % v1;
       norm.normalize();
       
       for (m = 0; m < 3; m++)
         {
         switch(m)
           {
           case 0:
             n = n0;
           break;
           case 1:
             n = n1;
           break;
           case 2:
             n = n2;
           break;
           default:
           break;
           }
         //now, we need to add nodes and elements at given heights
         if (stack[n]->max == 0)
           {
           if (shared.Is_In_List(n))
             {
             nstore[0] = nstore[1] = nstore[2] = 0.0;
             int tl = 0;
             //need to change norm to avg of facets listed
             for (k = 0; k < geom->ngb; k++)
               {
               l = tempv[n][k];

               if (l < 0)
                 continue; 

               v0 = Vector(geom->g_vert[geom->g_facet[l].nodes[0]],geom->g_vert[geom->g_facet[l].nodes[1]]);
               v1 = Vector(geom->g_vert[geom->g_facet[l].nodes[0]],geom->g_vert[geom->g_facet[l].nodes[2]]);

               normv = v0 % v1;
               normv.normalize();

               nstore[0] += normv[0];
               nstore[1] += normv[1];
               nstore[2] += normv[2];
               tl++;
               }

             //now, reset normt
             normt = Vector(nstore[0]/tl,nstore[1]/tl,nstore[2]/tl);
             normt.normalize();
             }
           else
             {
             normt = Vector(norm[0],norm[1],norm[2]);
             }
           pv = geom->g_vert[n];
           for (k = 0; k < viscous; k++)
             {
             normv = normt*(vspace*pow(factor,k));
             pv += Point(normv[0],normv[1],normv[2]);

             //add node
             if (!overset)
               {
               finalmesh->node[finalmesh->nn].vert = pv;
               stack[n]->Add_To_List(finalmesh->nn);
               if (k == viscous-1)
                 bdnodes.Add_To_List(finalmesh->nn);
               finalmesh->nn++;
               }
             else
               {
               ovst->node[ovst->nn].vert = pv;
               stack[n]->Add_To_List(ovst->nn);
               ovst->nn++;
               }
             }
           }
         }
       }
     }

   shared.destruct();
   for (i = 0; i < geom->n_gnodes; i++)
     delete [] tempv[i];
   delete [] tempv;

  if (!overset)
    {
    //set up map for bd nodes
    map = new int[finalmesh->nn+geom->n_gnodes];
    //be sure init to -1
    for (k = 0; k < (finalmesh->nn+geom->n_gnodes); k++)
      map[k] = -4;

    //put bd nodes in map
    for (k = 0; k < bdnodes.max; k++)
      map[bdnodes.list[k]] = -2;

    //add in geom facets to tri
    ngt = 0;

  
    fprintf(out_f,"\nAllocating for geom boundary elements.\n");
    fflush(out_f);
  
    if ((finalmesh->nelems+geom->n_gfacets) > finalmesh->element_dim)
      my_mem(finalmesh->element_dim, finalmesh->nelems, &(finalmesh->element), (finalmesh->nelems+geom->n_gfacets), ELEMENT_CHUNK);
  
    //since all elements are one face, set up nf and f_n now
    for (k = finalmesh->nelems; k < (finalmesh->nelems+geom->n_gfacets); k++)
      {
      finalmesh->element[k].Initialize();
      finalmesh->element[k].nf = 1;
      finalmesh->element[k].f_n = new List*[1];
      finalmesh->element[k].f_n[0] = new List();
      }
    
    fprintf(out_f,"\nDetermining geometry tolerance and storing facets as boundary objects.\n");
    fflush(out_f);
      
    //reset tol for MIN op
    gtol = 1.0e+20;

    nonviscfacets = 0;
  
    //loop through geom facets, noting all unique nodes at this point
    //set up geom as elements, then tris as elements
    for (i=1; i <= geom->ngb; i++)
      {
      for (j=geom->n_begin[i]; j < geom->n_end[i]; j++)
        {
        n0 = geom->g_facet[j].nodes[0];
        n1 = geom->g_facet[j].nodes[1];
        n2 = geom->g_facet[j].nodes[2];
      
        gtol = MIN(gtol,distance(geom->g_vert[n0], geom->g_vert[n1]));
        gtol = MIN(gtol,distance(geom->g_vert[n1], geom->g_vert[n2]));
        gtol = MIN(gtol,distance(geom->g_vert[n2], geom->g_vert[n0]));
      
        ngt++;
        //now, need to create tris
        /*if ((finalmesh->nelems+1) > finalmesh->element_dim)
          {
          my_mem(finalmesh->element_dim, finalmesh->nelems, &(finalmesh->element), (1+finalmesh->nelems), ELEMENT_CHUNK);
      
          //init elems
          for (n = finalmesh->nelems; n < (finalmesh->nelems+1); n++)
            {
            finalmesh->element[n].Initialize();
            finalmesh->element[n].nf = 1;
            finalmesh->element[n].f_n = new List*[1];
            finalmesh->element[n].f_n[0] = new List();
            }
          } */ 

        //take the resulting triangles and set up tri
        finalmesh->element[finalmesh->nelems].f_n[0]->Add_To_List(n0+finalmesh->nn);
        finalmesh->element[finalmesh->nelems].f_n[0]->Add_To_List(n1+finalmesh->nn);
        finalmesh->element[finalmesh->nelems].f_n[0]->Add_To_List(n2+finalmesh->nn);

        if (!nonviscbd.Is_In_List(i))
          {
          map[n0+finalmesh->nn] = -3;
          map[n1+finalmesh->nn] = -3;
          map[n2+finalmesh->nn] = -3;
          }
        else
          {
          map[n0+finalmesh->nn] = -1;
          map[n1+finalmesh->nn] = -1;
          map[n2+finalmesh->nn] = -1;
          nonviscfacets++;
          }
      
        finalmesh->boundary[i-1].elist->Add_To_List(finalmesh->nelems); 
        //inc nelems
        finalmesh->nelems++;
        }
      }

    fprintf(out_f,"\n#Elem = %d\n",finalmesh->nelems);
    fflush(out_f); 

    //finalize gtol, print out so can be sent into restart
    gtol *= 0.25;
  
    fprintf(out_f,"\nGeometry tolerance is %16.10e.  Please enter after restart.\n",gtol);
    fflush(out_f);
  
    //set up tol as 0.25 smn, print out so can be sent into restart
    tol = 0.25*smn;
  
    fprintf(out_f,"\nInternal boundary tolerance is %16.10e.  Please enter after restart.\n",tol);
    fflush(out_f);
    
    if (ngt != geom->n_gfacets)
      {
      fprintf(out_f,"\nNumber of geometry facets doesn't match number to be written to UCD file. Exiting....\n");
      fflush(out_f);
      exit(0);
      }
    }

   //int currne = finalmesh->nelems; //save for looping later
   int numfaces = 0;
   List directnode;
   directnode.construct();
   List indirectnode;
   indirectnode.construct();

   if (!overset)
     {
     my_mem(finalmesh->element_dim, finalmesh->nelems, &(finalmesh->element), finalmesh->nelems+((viscous+1)*geom->n_gfacets), ELEMENT_CHUNK);
     for (k = finalmesh->nelems; k < finalmesh->nelems+(viscous+1)*geom->n_gfacets; k++)
       finalmesh->element[k].Initialize();
     }

   if (overset && ovst->nelems+((viscous+1)*geom->n_gfacets) > ovst->element_dim)
     {
     my_mem(ovst->element_dim, ovst->nelems, &(ovst->element), ovst->nelems+((viscous+1)*geom->n_gfacets), ELEMENT_CHUNK);
     for (k = ovst->nelems; k < ovst->nelems+(viscous+1)*geom->n_gfacets; k++)
       ovst->element[k].Initialize();
     }

   //in order to assure that we have real elem created before geom/stitch/visc stitch, will create prisms, then bd faces
   //now, that we have created stacks, we need to create elements, as well as add to stitching bd
   for (i = 0; i < nvbd; i++)
     {
     for (j=geom->n_begin[vbd[i]]; j < geom->n_end[vbd[i]]; j++)
       {
       n0 = geom->g_facet[j].nodes[0];
       n1 = geom->g_facet[j].nodes[1];
       n2 = geom->g_facet[j].nodes[2];

       if (keepfacet.Is_In_List(j))
       {
       directnode.Redimension(0);
       indirectnode.Redimension(0);
       if (stack[n0]->max == 0 && stack[n1]->max == 0 && stack[n2]->max == 0)
         {
         fprintf(stderr,"\nYou have a facet %d adjacent to a viscous bd that has no viscous nabor (%lf, %lf, %lf).  Exiting....\n",j,geom->g_vert[n0][0],geom->g_vert[n0][1],geom->g_vert[n0][2]);
         fflush(stderr);
         exit(0);
         }
       else if (stack[n0]->max == 0 && stack[n1]->max ==0)
         {
         numfaces = 2*viscous+2;
         directnode.Add_To_List(n1);
         directnode.Add_To_List(n0);
         indirectnode.Add_To_List(n2);
         }
       else if (stack[n0]->max == 0 && stack[n2]->max ==0)
         {
         numfaces = 2*viscous+2;
         directnode.Add_To_List(n0);
         directnode.Add_To_List(n2);
         indirectnode.Add_To_List(n1);
         }
       else if (stack[n1]->max == 0 && stack[n2]->max ==0)
         {
         numfaces = 2*viscous+2;
         directnode.Add_To_List(n2);
         directnode.Add_To_List(n1);
         indirectnode.Add_To_List(n0);
         }
       else if (stack[n0]->max == 0)
         {
         numfaces = 3*viscous+2;
         directnode.Add_To_List(n0);
         indirectnode.Add_To_List(n2);
         indirectnode.Add_To_List(n1);
         }
       else if (stack[n1]->max == 0)
         {
         numfaces = 3*viscous+2;
         directnode.Add_To_List(n1);
         indirectnode.Add_To_List(n0);
         indirectnode.Add_To_List(n2);
         }
       else if (stack[n2]->max == 0)
         {
         numfaces = 3*viscous+2;
         directnode.Add_To_List(n2);
         indirectnode.Add_To_List(n1);
         indirectnode.Add_To_List(n0);
         }
       if (!overset)
       {
       finalmesh->element[finalmesh->nelems].nf = numfaces;
       finalmesh->element[finalmesh->nelems].f_n = new List*[numfaces];
       for (l = 0; l < numfaces; l++)
         finalmesh->element[finalmesh->nelems].f_n[l] = new List();
       //use stacks to create elements
       if (directnode.max == 1)
         {
         finalmesh->element[finalmesh->nelems].f_n[0]->Add_To_List(directnode.list[0]+finalmesh->nn);
         finalmesh->element[finalmesh->nelems].f_n[0]->Add_To_List(indirectnode.list[0]+finalmesh->nn);
         finalmesh->element[finalmesh->nelems].f_n[0]->Add_To_List(indirectnode.list[1]+finalmesh->nn);

         finalmesh->element[finalmesh->nelems].f_n[1]->Add_To_List(directnode.list[0]+finalmesh->nn);
         finalmesh->element[finalmesh->nelems].f_n[1]->Add_To_List(stack[indirectnode.list[1]]->list[stack[indirectnode.list[1]]->max-1]);
         finalmesh->element[finalmesh->nelems].f_n[1]->Add_To_List(stack[indirectnode.list[0]]->list[stack[indirectnode.list[0]]->max-1]);

         finalmesh->element[finalmesh->nelems].f_n[2]->Add_To_List(directnode.list[0]+finalmesh->nn);
         finalmesh->element[finalmesh->nelems].f_n[2]->Add_To_List(stack[indirectnode.list[0]]->list[0]);
         finalmesh->element[finalmesh->nelems].f_n[2]->Add_To_List(indirectnode.list[0]+finalmesh->nn);

         finalmesh->element[finalmesh->nelems].f_n[3]->Add_To_List(directnode.list[0]+finalmesh->nn);
         finalmesh->element[finalmesh->nelems].f_n[3]->Add_To_List(indirectnode.list[1]+finalmesh->nn);
         finalmesh->element[finalmesh->nelems].f_n[3]->Add_To_List(stack[indirectnode.list[1]]->list[0]);

         finalmesh->element[finalmesh->nelems].f_n[4]->Add_To_List(indirectnode.list[1]+finalmesh->nn);
         finalmesh->element[finalmesh->nelems].f_n[4]->Add_To_List(indirectnode.list[0]+finalmesh->nn);
         finalmesh->element[finalmesh->nelems].f_n[4]->Add_To_List(stack[indirectnode.list[0]]->list[0]);
         finalmesh->element[finalmesh->nelems].f_n[4]->Add_To_List(stack[indirectnode.list[1]]->list[0]);

         l = 5; //counter
         for (k = 0; k < viscous-1; k++)
           {
           finalmesh->element[finalmesh->nelems].f_n[l]->Add_To_List(directnode.list[0]+finalmesh->nn);
           finalmesh->element[finalmesh->nelems].f_n[l]->Add_To_List(stack[indirectnode.list[0]]->list[k+1]);
           finalmesh->element[finalmesh->nelems].f_n[l]->Add_To_List(stack[indirectnode.list[0]]->list[k]);
           l++;

           finalmesh->element[finalmesh->nelems].f_n[l]->Add_To_List(directnode.list[0]+finalmesh->nn);
           finalmesh->element[finalmesh->nelems].f_n[l]->Add_To_List(stack[indirectnode.list[1]]->list[k]);
           finalmesh->element[finalmesh->nelems].f_n[l]->Add_To_List(stack[indirectnode.list[1]]->list[k+1]);
           l++;

           finalmesh->element[finalmesh->nelems].f_n[l]->Add_To_List(stack[indirectnode.list[1]]->list[k]);
           finalmesh->element[finalmesh->nelems].f_n[l]->Add_To_List(stack[indirectnode.list[0]]->list[k]);
           finalmesh->element[finalmesh->nelems].f_n[l]->Add_To_List(stack[indirectnode.list[0]]->list[k+1]);
           finalmesh->element[finalmesh->nelems].f_n[l]->Add_To_List(stack[indirectnode.list[1]]->list[k+1]);
           l++;
           }
         finalmesh->nelems++;
         }
       if (directnode.max == 2)
         {
         finalmesh->element[finalmesh->nelems].f_n[0]->Add_To_List(directnode.list[0]+finalmesh->nn);
         finalmesh->element[finalmesh->nelems].f_n[0]->Add_To_List(directnode.list[1]+finalmesh->nn);
         finalmesh->element[finalmesh->nelems].f_n[0]->Add_To_List(indirectnode.list[0]+finalmesh->nn);

         finalmesh->element[finalmesh->nelems].f_n[1]->Add_To_List(directnode.list[0]+finalmesh->nn);
         finalmesh->element[finalmesh->nelems].f_n[1]->Add_To_List(stack[indirectnode.list[0]]->list[stack[indirectnode.list[0]]->max-1]);
         finalmesh->element[finalmesh->nelems].f_n[1]->Add_To_List(directnode.list[1]+finalmesh->nn);

         finalmesh->element[finalmesh->nelems].f_n[2]->Add_To_List(directnode.list[0]+finalmesh->nn);
         finalmesh->element[finalmesh->nelems].f_n[2]->Add_To_List(indirectnode.list[0]+finalmesh->nn);
         finalmesh->element[finalmesh->nelems].f_n[2]->Add_To_List(stack[indirectnode.list[0]]->list[0]);

         finalmesh->element[finalmesh->nelems].f_n[3]->Add_To_List(stack[indirectnode.list[0]]->list[0]);
         finalmesh->element[finalmesh->nelems].f_n[3]->Add_To_List(indirectnode.list[0]+finalmesh->nn);
         finalmesh->element[finalmesh->nelems].f_n[3]->Add_To_List(directnode.list[1]+finalmesh->nn);

         l = 4; //counter
         for (k = 0; k < viscous-1; k++)
           {
           finalmesh->element[finalmesh->nelems].f_n[l]->Add_To_List(directnode.list[0]+finalmesh->nn);
           finalmesh->element[finalmesh->nelems].f_n[l]->Add_To_List(stack[indirectnode.list[0]]->list[k]);
           finalmesh->element[finalmesh->nelems].f_n[l]->Add_To_List(stack[indirectnode.list[0]]->list[k+1]);
           l++;

           finalmesh->element[finalmesh->nelems].f_n[l]->Add_To_List(stack[indirectnode.list[0]]->list[k+1]);
           finalmesh->element[finalmesh->nelems].f_n[l]->Add_To_List(stack[indirectnode.list[0]]->list[k]);
           finalmesh->element[finalmesh->nelems].f_n[l]->Add_To_List(directnode.list[1]+finalmesh->nn);
           l++;
           }
         finalmesh->nelems++;
         }
       }
       if (overset)
       {
       ovst->element[ovst->nelems].nf = numfaces;
       ovst->element[ovst->nelems].f_n = new List*[numfaces];
       for (l = 0; l < numfaces; l++)
         ovst->element[ovst->nelems].f_n[l] = new List();
       //use stacks to create elements
       if (directnode.max == 1)
         {
         ovst->element[ovst->nelems].f_n[0]->Add_To_List(omap[directnode.list[0]]);
         ovst->element[ovst->nelems].f_n[0]->Add_To_List(omap[indirectnode.list[0]]);
         ovst->element[ovst->nelems].f_n[0]->Add_To_List(omap[indirectnode.list[1]]);

         ovst->element[ovst->nelems].f_n[1]->Add_To_List(omap[directnode.list[0]]);
         ovst->element[ovst->nelems].f_n[1]->Add_To_List(stack[indirectnode.list[1]]->list[stack[indirectnode.list[1]]->max-1]);
         ovst->element[ovst->nelems].f_n[1]->Add_To_List(stack[indirectnode.list[0]]->list[stack[indirectnode.list[0]]->max-1]);

         ovst->element[ovst->nelems].f_n[2]->Add_To_List(omap[directnode.list[0]]);
         ovst->element[ovst->nelems].f_n[2]->Add_To_List(stack[indirectnode.list[0]]->list[0]);
         ovst->element[ovst->nelems].f_n[2]->Add_To_List(omap[indirectnode.list[0]]);

         ovst->element[ovst->nelems].f_n[3]->Add_To_List(omap[directnode.list[0]]);
         ovst->element[ovst->nelems].f_n[3]->Add_To_List(omap[indirectnode.list[1]]);
         ovst->element[ovst->nelems].f_n[3]->Add_To_List(stack[indirectnode.list[1]]->list[0]);

         ovst->element[ovst->nelems].f_n[4]->Add_To_List(omap[indirectnode.list[1]]);
         ovst->element[ovst->nelems].f_n[4]->Add_To_List(omap[indirectnode.list[0]]);
         ovst->element[ovst->nelems].f_n[4]->Add_To_List(stack[indirectnode.list[0]]->list[0]);
         ovst->element[ovst->nelems].f_n[4]->Add_To_List(stack[indirectnode.list[1]]->list[0]);

         l = 5; //counter
         for (k = 0; k < viscous-1; k++)
           {
           ovst->element[ovst->nelems].f_n[l]->Add_To_List(omap[directnode.list[0]]);
           ovst->element[ovst->nelems].f_n[l]->Add_To_List(stack[indirectnode.list[0]]->list[k+1]);
           ovst->element[ovst->nelems].f_n[l]->Add_To_List(stack[indirectnode.list[0]]->list[k]);
           l++;

           ovst->element[ovst->nelems].f_n[l]->Add_To_List(omap[directnode.list[0]]);
           ovst->element[ovst->nelems].f_n[l]->Add_To_List(stack[indirectnode.list[1]]->list[k]);
           ovst->element[ovst->nelems].f_n[l]->Add_To_List(stack[indirectnode.list[1]]->list[k+1]);
           l++;

           ovst->element[ovst->nelems].f_n[l]->Add_To_List(stack[indirectnode.list[1]]->list[k]);
           ovst->element[ovst->nelems].f_n[l]->Add_To_List(stack[indirectnode.list[0]]->list[k]);
           ovst->element[ovst->nelems].f_n[l]->Add_To_List(stack[indirectnode.list[0]]->list[k+1]);
           ovst->element[ovst->nelems].f_n[l]->Add_To_List(stack[indirectnode.list[1]]->list[k+1]);
           l++;
           }
         ovst->nelems++;
         }
       if (directnode.max == 2)
         {
         ovst->element[ovst->nelems].f_n[0]->Add_To_List(omap[directnode.list[0]]);
         ovst->element[ovst->nelems].f_n[0]->Add_To_List(omap[directnode.list[1]]);
         ovst->element[ovst->nelems].f_n[0]->Add_To_List(omap[indirectnode.list[0]]);

         ovst->element[ovst->nelems].f_n[1]->Add_To_List(omap[directnode.list[0]]);
         ovst->element[ovst->nelems].f_n[1]->Add_To_List(stack[indirectnode.list[0]]->list[stack[indirectnode.list[0]]->max-1]);
         ovst->element[ovst->nelems].f_n[1]->Add_To_List(omap[directnode.list[1]]);

         ovst->element[ovst->nelems].f_n[2]->Add_To_List(omap[directnode.list[0]]);
         ovst->element[ovst->nelems].f_n[2]->Add_To_List(omap[indirectnode.list[0]]);
         ovst->element[ovst->nelems].f_n[2]->Add_To_List(stack[indirectnode.list[0]]->list[0]);

         ovst->element[ovst->nelems].f_n[3]->Add_To_List(stack[indirectnode.list[0]]->list[0]);
         ovst->element[ovst->nelems].f_n[3]->Add_To_List(omap[indirectnode.list[0]]);
         ovst->element[ovst->nelems].f_n[3]->Add_To_List(omap[directnode.list[1]]);

         l = 4; //counter
         for (k = 0; k < viscous-1; k++)
           {
           ovst->element[ovst->nelems].f_n[l]->Add_To_List(omap[directnode.list[0]]);
           ovst->element[ovst->nelems].f_n[l]->Add_To_List(stack[indirectnode.list[0]]->list[k]);
           ovst->element[ovst->nelems].f_n[l]->Add_To_List(stack[indirectnode.list[0]]->list[k+1]);
           l++;

           ovst->element[ovst->nelems].f_n[l]->Add_To_List(stack[indirectnode.list[0]]->list[k+1]);
           ovst->element[ovst->nelems].f_n[l]->Add_To_List(stack[indirectnode.list[0]]->list[k]);
           ovst->element[ovst->nelems].f_n[l]->Add_To_List(omap[directnode.list[1]]);
           l++;
           }
         ovst->nelems++;
         }
       }
       }
       else
       {
       if (!overset)
       {
       //create base element
       finalmesh->element[finalmesh->nelems].nf = 5;
       finalmesh->element[finalmesh->nelems].f_n = new List*[5];
       for (l = 0; l < 5; l++)
         finalmesh->element[finalmesh->nelems].f_n[l] = new List();

       //add in faces in CGNS order
       finalmesh->element[finalmesh->nelems].f_n[0]->Add_To_List(n0+finalmesh->nn);
       finalmesh->element[finalmesh->nelems].f_n[0]->Add_To_List(n1+finalmesh->nn);
       finalmesh->element[finalmesh->nelems].f_n[0]->Add_To_List(stack[n1]->list[0]);
       finalmesh->element[finalmesh->nelems].f_n[0]->Add_To_List(stack[n0]->list[0]);

       finalmesh->element[finalmesh->nelems].f_n[1]->Add_To_List(n1+finalmesh->nn);
       finalmesh->element[finalmesh->nelems].f_n[1]->Add_To_List(n2+finalmesh->nn);
       finalmesh->element[finalmesh->nelems].f_n[1]->Add_To_List(stack[n2]->list[0]);
       finalmesh->element[finalmesh->nelems].f_n[1]->Add_To_List(stack[n1]->list[0]);

       finalmesh->element[finalmesh->nelems].f_n[2]->Add_To_List(n2+finalmesh->nn);
       finalmesh->element[finalmesh->nelems].f_n[2]->Add_To_List(n0+finalmesh->nn);
       finalmesh->element[finalmesh->nelems].f_n[2]->Add_To_List(stack[n0]->list[0]);
       finalmesh->element[finalmesh->nelems].f_n[2]->Add_To_List(stack[n2]->list[0]);

       finalmesh->element[finalmesh->nelems].f_n[3]->Add_To_List(n0+finalmesh->nn);
       finalmesh->element[finalmesh->nelems].f_n[3]->Add_To_List(n2+finalmesh->nn);
       finalmesh->element[finalmesh->nelems].f_n[3]->Add_To_List(n1+finalmesh->nn);

       finalmesh->element[finalmesh->nelems].f_n[4]->Add_To_List(stack[n0]->list[0]);
       finalmesh->element[finalmesh->nelems].f_n[4]->Add_To_List(stack[n1]->list[0]);
       finalmesh->element[finalmesh->nelems].f_n[4]->Add_To_List(stack[n2]->list[0]);

       finalmesh->nelems++;
       //use stacks to create rest of elements
       for (k = 0; k < viscous-1; k++)
         {
         finalmesh->element[finalmesh->nelems].nf = 5;
         finalmesh->element[finalmesh->nelems].f_n = new List*[5];
         for (l = 0; l < 5; l++)
           finalmesh->element[finalmesh->nelems].f_n[l] = new List();

         //add in faces in CGNS order
         finalmesh->element[finalmesh->nelems].f_n[0]->Add_To_List(stack[n0]->list[k]);
         finalmesh->element[finalmesh->nelems].f_n[0]->Add_To_List(stack[n1]->list[k]);
         finalmesh->element[finalmesh->nelems].f_n[0]->Add_To_List(stack[n1]->list[k+1]);
         finalmesh->element[finalmesh->nelems].f_n[0]->Add_To_List(stack[n0]->list[k+1]);

         finalmesh->element[finalmesh->nelems].f_n[1]->Add_To_List(stack[n1]->list[k]);
         finalmesh->element[finalmesh->nelems].f_n[1]->Add_To_List(stack[n2]->list[k]);
         finalmesh->element[finalmesh->nelems].f_n[1]->Add_To_List(stack[n2]->list[k+1]);
         finalmesh->element[finalmesh->nelems].f_n[1]->Add_To_List(stack[n1]->list[k+1]);

         finalmesh->element[finalmesh->nelems].f_n[2]->Add_To_List(stack[n2]->list[k]);
         finalmesh->element[finalmesh->nelems].f_n[2]->Add_To_List(stack[n0]->list[k]);
         finalmesh->element[finalmesh->nelems].f_n[2]->Add_To_List(stack[n0]->list[k+1]);
         finalmesh->element[finalmesh->nelems].f_n[2]->Add_To_List(stack[n2]->list[k+1]);

         finalmesh->element[finalmesh->nelems].f_n[3]->Add_To_List(stack[n0]->list[k]);
         finalmesh->element[finalmesh->nelems].f_n[3]->Add_To_List(stack[n2]->list[k]);
         finalmesh->element[finalmesh->nelems].f_n[3]->Add_To_List(stack[n1]->list[k]);

         finalmesh->element[finalmesh->nelems].f_n[4]->Add_To_List(stack[n0]->list[k+1]);
         finalmesh->element[finalmesh->nelems].f_n[4]->Add_To_List(stack[n1]->list[k+1]);
         finalmesh->element[finalmesh->nelems].f_n[4]->Add_To_List(stack[n2]->list[k+1]);

         finalmesh->nelems++;
         }
       }
       if (overset)
       {
       //create base element
       ovst->element[ovst->nelems].nf = 5;
       ovst->element[ovst->nelems].f_n = new List*[5];
       for (l = 0; l < 5; l++)
         ovst->element[ovst->nelems].f_n[l] = new List();

       //add in faces in CGNS order
       ovst->element[ovst->nelems].f_n[0]->Add_To_List(omap[n0]);
       ovst->element[ovst->nelems].f_n[0]->Add_To_List(omap[n1]);
       ovst->element[ovst->nelems].f_n[0]->Add_To_List(stack[n1]->list[0]);
       ovst->element[ovst->nelems].f_n[0]->Add_To_List(stack[n0]->list[0]);

       ovst->element[ovst->nelems].f_n[1]->Add_To_List(omap[n1]);
       ovst->element[ovst->nelems].f_n[1]->Add_To_List(omap[n2]);
       ovst->element[ovst->nelems].f_n[1]->Add_To_List(stack[n2]->list[0]);
       ovst->element[ovst->nelems].f_n[1]->Add_To_List(stack[n1]->list[0]);

       ovst->element[ovst->nelems].f_n[2]->Add_To_List(omap[n2]);
       ovst->element[ovst->nelems].f_n[2]->Add_To_List(omap[n0]);
       ovst->element[ovst->nelems].f_n[2]->Add_To_List(stack[n0]->list[0]);
       ovst->element[ovst->nelems].f_n[2]->Add_To_List(stack[n2]->list[0]);

       ovst->element[ovst->nelems].f_n[3]->Add_To_List(omap[n0]);
       ovst->element[ovst->nelems].f_n[3]->Add_To_List(omap[n2]);
       ovst->element[ovst->nelems].f_n[3]->Add_To_List(omap[n1]);

       ovst->element[ovst->nelems].f_n[4]->Add_To_List(stack[n0]->list[0]);
       ovst->element[ovst->nelems].f_n[4]->Add_To_List(stack[n1]->list[0]);
       ovst->element[ovst->nelems].f_n[4]->Add_To_List(stack[n2]->list[0]);

       ovst->nelems++;
       //use stacks to create rest of elements
       for (k = 0; k < viscous-1; k++)
         {
         ovst->element[ovst->nelems].nf = 5;
         ovst->element[ovst->nelems].f_n = new List*[5];
         for (l = 0; l < 5; l++)
           ovst->element[ovst->nelems].f_n[l] = new List();

         //add in faces in CGNS order
         ovst->element[ovst->nelems].f_n[0]->Add_To_List(stack[n0]->list[k]);
         ovst->element[ovst->nelems].f_n[0]->Add_To_List(stack[n1]->list[k]);
         ovst->element[ovst->nelems].f_n[0]->Add_To_List(stack[n1]->list[k+1]);
         ovst->element[ovst->nelems].f_n[0]->Add_To_List(stack[n0]->list[k+1]);

         ovst->element[ovst->nelems].f_n[1]->Add_To_List(stack[n1]->list[k]);
         ovst->element[ovst->nelems].f_n[1]->Add_To_List(stack[n2]->list[k]);
         ovst->element[ovst->nelems].f_n[1]->Add_To_List(stack[n2]->list[k+1]);
         ovst->element[ovst->nelems].f_n[1]->Add_To_List(stack[n1]->list[k+1]);

         ovst->element[ovst->nelems].f_n[2]->Add_To_List(stack[n2]->list[k]);
         ovst->element[ovst->nelems].f_n[2]->Add_To_List(stack[n0]->list[k]);
         ovst->element[ovst->nelems].f_n[2]->Add_To_List(stack[n0]->list[k+1]);
         ovst->element[ovst->nelems].f_n[2]->Add_To_List(stack[n2]->list[k+1]);

         ovst->element[ovst->nelems].f_n[3]->Add_To_List(stack[n0]->list[k]);
         ovst->element[ovst->nelems].f_n[3]->Add_To_List(stack[n2]->list[k]);
         ovst->element[ovst->nelems].f_n[3]->Add_To_List(stack[n1]->list[k]);

         ovst->element[ovst->nelems].f_n[4]->Add_To_List(stack[n0]->list[k+1]);
         ovst->element[ovst->nelems].f_n[4]->Add_To_List(stack[n1]->list[k+1]);
         ovst->element[ovst->nelems].f_n[4]->Add_To_List(stack[n2]->list[k+1]);

         ovst->nelems++;
         }
       }
       }
       }
     }
  //just doing bd faces
  for (i = 0; i < nvbd; i++)
     {
     for (j=geom->n_begin[vbd[i]]; j < geom->n_end[vbd[i]]; j++)
       {
       n0 = geom->g_facet[j].nodes[0];
       n1 = geom->g_facet[j].nodes[1];
       n2 = geom->g_facet[j].nodes[2];

       if (keepfacet.Is_In_List(j))
       {
       directnode.Redimension(0);
       indirectnode.Redimension(0);
       if (stack[n0]->max == 0 && stack[n1]->max == 0 && stack[n2]->max == 0)
         {
         fprintf(stderr,"\nYou have a facet %d adjacent to a viscous bd that has no viscous nabor (%lf, %lf, %lf).  Exiting....\n",j,geom->g_vert[n0][0],geom->g_vert[n0][1],geom->g_vert[n0][2]);
         fflush(stderr);
         exit(0);
         }
       else if (stack[n0]->max == 0 && stack[n1]->max ==0)
         {
         numfaces = 2*viscous+2;
         directnode.Add_To_List(n1);
         directnode.Add_To_List(n0);
         indirectnode.Add_To_List(n2);
         }
       else if (stack[n0]->max == 0 && stack[n2]->max ==0)
         {
         numfaces = 2*viscous+2;
         directnode.Add_To_List(n0);
         directnode.Add_To_List(n2);
         indirectnode.Add_To_List(n1);
         }
       else if (stack[n1]->max == 0 && stack[n2]->max ==0)
         {
         numfaces = 2*viscous+2;
         directnode.Add_To_List(n2);
         directnode.Add_To_List(n1);
         indirectnode.Add_To_List(n0);
         }
       else if (stack[n0]->max == 0)
         {
         numfaces = 3*viscous+2;
         directnode.Add_To_List(n0);
         indirectnode.Add_To_List(n2);
         indirectnode.Add_To_List(n1);
         }
       else if (stack[n1]->max == 0)
         {
         numfaces = 3*viscous+2;
         directnode.Add_To_List(n1);
         indirectnode.Add_To_List(n0);
         indirectnode.Add_To_List(n2);
         }
       else if (stack[n2]->max == 0)
         {
         numfaces = 3*viscous+2;
         directnode.Add_To_List(n2);
         indirectnode.Add_To_List(n1);
         indirectnode.Add_To_List(n0);
         }
       //use stacks to create elements
       if (!overset)
       {
       if (directnode.max == 1)
         {
         //finally, put top face on stitching bd
         finalmesh->element[finalmesh->nelems].nf = 1;
         finalmesh->element[finalmesh->nelems].f_n = new List*[1];
         finalmesh->element[finalmesh->nelems].f_n[0] = new List();
         finalmesh->element[finalmesh->nelems].f_n[0]->Add_To_List(directnode.list[0]+finalmesh->nn);
         finalmesh->element[finalmesh->nelems].f_n[0]->Add_To_List(stack[indirectnode.list[1]]->list[stack[indirectnode.list[1]]->max-1]);
         finalmesh->element[finalmesh->nelems].f_n[0]->Add_To_List(stack[indirectnode.list[0]]->list[stack[indirectnode.list[0]]->max-1]);
         finalmesh->boundary[geom->ngb+1].elist->Add_To_List(finalmesh->nelems);
         tempnt++;
         finalmesh->nelems++;
         }
       if (directnode.max == 2)
         {
         //finally, put top face on stitching bd
         finalmesh->element[finalmesh->nelems].nf = 1;
         finalmesh->element[finalmesh->nelems].f_n = new List*[1];
         finalmesh->element[finalmesh->nelems].f_n[0] = new List();
         finalmesh->element[finalmesh->nelems].f_n[0]->Add_To_List(directnode.list[0]+finalmesh->nn);
         finalmesh->element[finalmesh->nelems].f_n[0]->Add_To_List(stack[indirectnode.list[0]]->list[stack[indirectnode.list[0]]->max-1]);
         finalmesh->element[finalmesh->nelems].f_n[0]->Add_To_List(directnode.list[1]+finalmesh->nn);
         finalmesh->boundary[geom->ngb+1].elist->Add_To_List(finalmesh->nelems);
         tempnt++;
         finalmesh->nelems++;
         }
       }
       if (overset)
       {
       if (directnode.max == 1)
         {
         //finally, put top face on stitching bd
         ovst->element[ovst->nelems].nf = 1;
         ovst->element[ovst->nelems].f_n = new List*[1];
         ovst->element[ovst->nelems].f_n[0] = new List();
         ovst->element[ovst->nelems].f_n[0]->Add_To_List(omap[directnode.list[0]]);
         ovst->element[ovst->nelems].f_n[0]->Add_To_List(stack[indirectnode.list[1]]->list[stack[indirectnode.list[1]]->max-1]);
         ovst->element[ovst->nelems].f_n[0]->Add_To_List(stack[indirectnode.list[0]]->list[stack[indirectnode.list[0]]->max-1]);
         ovst->boundary[nvbd].elist->Add_To_List(ovst->nelems);
         tempnt++;
         ovst->nelems++;
         }
       if (directnode.max == 2)
         {
         //finally, put top face on stitching bd
         ovst->element[ovst->nelems].nf = 1;
         ovst->element[ovst->nelems].f_n = new List*[1];
         ovst->element[ovst->nelems].f_n[0] = new List();
         ovst->element[ovst->nelems].f_n[0]->Add_To_List(omap[directnode.list[0]]);
         ovst->element[ovst->nelems].f_n[0]->Add_To_List(stack[indirectnode.list[0]]->list[stack[indirectnode.list[0]]->max-1]);
         ovst->element[ovst->nelems].f_n[0]->Add_To_List(omap[directnode.list[1]]);
         ovst->boundary[nvbd].elist->Add_To_List(ovst->nelems);
         tempnt++;
         ovst->nelems++;
         }
       }
       }
       else
       {
       if (!overset)
       {
       //finally, put top face on stitching bd
       finalmesh->element[finalmesh->nelems].nf = 1;
       finalmesh->element[finalmesh->nelems].f_n = new List*[1];
       finalmesh->element[finalmesh->nelems].f_n[0] = new List();
       finalmesh->element[finalmesh->nelems].f_n[0]->Add_To_List(stack[n0]->list[stack[n0]->max - 1]);
       finalmesh->element[finalmesh->nelems].f_n[0]->Add_To_List(stack[n1]->list[stack[n1]->max - 1]);
       finalmesh->element[finalmesh->nelems].f_n[0]->Add_To_List(stack[n2]->list[stack[n2]->max - 1]);
       finalmesh->boundary[geom->ngb+1].elist->Add_To_List(finalmesh->nelems);
       tempnt++;
       finalmesh->nelems++;
       }
       if (overset)
       {
       //finally, put top face on stitching bd
       ovst->element[ovst->nelems].nf = 1;
       ovst->element[ovst->nelems].f_n = new List*[1];
       ovst->element[ovst->nelems].f_n[0] = new List();
       ovst->element[ovst->nelems].f_n[0]->Add_To_List(stack[n0]->list[stack[n0]->max - 1]);
       ovst->element[ovst->nelems].f_n[0]->Add_To_List(stack[n1]->list[stack[n1]->max - 1]);
       ovst->element[ovst->nelems].f_n[0]->Add_To_List(stack[n2]->list[stack[n2]->max - 1]);
       ovst->boundary[nvbd].elist->Add_To_List(ovst->nelems);
       tempnt++;
       ovst->nelems++;
       }
       }
       }
     }

  //add in geom nodes so map is right
  for (i = 0; i < geom->n_gnodes && !overset; i++)
    {
    if ((finalmesh->nn + 1) > finalmesh->node_dim)
      my_mem(finalmesh->node_dim, finalmesh->nn, &(finalmesh->node), (finalmesh->nn + 1), NODE_CHUNK);
      
    finalmesh->node[finalmesh->nn].vert = geom->g_vert[i];
    
    //reset nn
    finalmesh->nn++;
    }

  fprintf(out_f,"\n#Elem = %d\n",finalmesh->nelems);
  fflush(out_f);

  bdnodes.destruct();
  directnode.destruct();
  indirectnode.destruct();
  delete [] omap;

  //no need to do this as edges of viscous bd now eased off via slant elements
  /*//make octree-sorted element list of these
  Octree_Storage *Octree_root1;
  POLY_ELEMENT *pntr1 = 0;
  PList *elist1 = 0;
  double stol1[3]; //will hold slop factor always
  double ptt[3]; //will hold pt for octree retrieve
  
  //reset points
  p0 = Point(1.0e20, 1.0e20, 1.0e20);
  p1 = Point(-1.0e20, -1.0e20, -1.0e20);
  
  int counter = 0;
  //find extents for viscous layers
  for (i = currne; i < finalmesh->nelems; i++)
     {
     if (finalmesh->element[i].nf == 1)
       continue;
     //count number elements
     counter++;
     for (j = 0; j < finalmesh->element[i].nf; j++)
       {
       for (k = 0; k < finalmesh->element[i].f_n[j]->max; k++)
         {
         p0[0] = MIN(p0[0],finalmesh->node[finalmesh->element[i].f_n[j]->list[k]].vert[0]);
         p0[1] = MIN(p0[1],finalmesh->node[finalmesh->element[i].f_n[j]->list[k]].vert[1]);
         p0[2] = MIN(p0[2],finalmesh->node[finalmesh->element[i].f_n[j]->list[k]].vert[2]);
         p1[0] = MAX(p1[0],finalmesh->node[finalmesh->element[i].f_n[j]->list[k]].vert[0]);
         p1[1] = MAX(p1[1],finalmesh->node[finalmesh->element[i].f_n[j]->list[k]].vert[1]);
         p1[2] = MAX(p1[2],finalmesh->node[finalmesh->element[i].f_n[j]->list[k]].vert[2]);
         }
       }
    }
    
  for (j = 0; j < 3; j++)
    {
    olo[j] = p0[j];
    ohi[j] = p1[j];
    }

  //set for mom
  Octree_Storage *dummy1=0;
  //set up tlo, thi
  double (*tlo1)[3], (*thi1)[3];
  //set to reuse and pass in
  void* *ptr1;
  //create storage
  Octree_root1 = new Octree_Storage((Octree_Storage*)dummy1,olo,ohi);
  
  //set up ptrs
  ptr1 = new void*[counter];
  tlo1 = new double[counter][3];
  thi1 = new double[counter][3];
  
  int counter1 = 0;
  //find extents for viscous layers
  for (i = currne; i < finalmesh->nelems; i++)
     {
     if (finalmesh->element[i].nf == 1)
       continue;
     //reset points
     p0 = Point(1.0e20, 1.0e20, 1.0e20);
     p1 = Point(-1.0e20, -1.0e20, -1.0e20);
     for (j = 0; j < finalmesh->element[i].nf; j++)
       {
       for (k = 0; k < finalmesh->element[i].f_n[j]->max; k++)
         {
         p0[0] = MIN(p0[0],finalmesh->node[finalmesh->element[i].f_n[j]->list[k]].vert[0]);
         p0[1] = MIN(p0[1],finalmesh->node[finalmesh->element[i].f_n[j]->list[k]].vert[1]);
         p0[2] = MIN(p0[2],finalmesh->node[finalmesh->element[i].f_n[j]->list[k]].vert[2]);
         p1[0] = MAX(p1[0],finalmesh->node[finalmesh->element[i].f_n[j]->list[k]].vert[0]);
         p1[1] = MAX(p1[1],finalmesh->node[finalmesh->element[i].f_n[j]->list[k]].vert[1]);
         p1[2] = MAX(p1[2],finalmesh->node[finalmesh->element[i].f_n[j]->list[k]].vert[2]);
         }
       }

    tlo1[counter1][0] = p0[0];
    tlo1[counter1][1] = p0[1];
    tlo1[counter1][2] = p0[2];
    thi1[counter1][0] = p1[0];
    thi1[counter1][1] = p1[1];
    thi1[counter1][2] = p1[2];
    ptr1[counter1] = (void*)&(finalmesh->element[i]);
    //count number elements
    counter1++;
    }
  
  //actually store facets
  Octree_root1->Store_In_Octree(counter1,ptr1,tlo1,thi1);

  //free mem
  delete[] ptr1;
  delete[] tlo1;
  delete[] thi1;
  
  //be sure pntr reset
  pntr1 = 0;
  
  //create new Plist
  elist1 = new PList();
  //dimension a priori, reset max within
  elist1->Redimension(counter);
  
  //set tolarance as lower of geom and inner tol
  tol = 0.1*fl;
 
   //finally, loop thru recently created to determine if any "sides" need to be on stitching bd (check for face 0 only and ignore)
   for (i = currne; i < finalmesh->nelems; i++)
     {
     if (finalmesh->element[i].nf == 1)
       continue;

     //find centroid
     nstore[0] = nstore[1] = nstore[2] = 0.0;
     for (j = 0; j < finalmesh->element[i].nf; j++)
       {
       for (k = 0; k < finalmesh->element[i].f_n[j]->max; k++)
         {
         nstore[0] += finalmesh->node[finalmesh->element[i].f_n[j]->list[k]].vert[0];
         nstore[1] += finalmesh->node[finalmesh->element[i].f_n[j]->list[k]].vert[1];
         nstore[2] += finalmesh->node[finalmesh->element[i].f_n[j]->list[k]].vert[2];
         }
       }

     nstore[0] /= 18.0;
     nstore[1] /= 18.0;
     nstore[2] /= 18.0;

     cent = Point(nstore[0],nstore[1],nstore[2]);

     //set level to -1 so returns all (only reason to go thru levels is to change slop factor) 
     l = -1;
     //now, redim elist
     elist1->max = 0;
     nstore[0] = nstore[1] = nstore[2] = 0.0;
     for (j = 0; j < finalmesh->element[i].nf; j++)
       {
       if (finalmesh->element[i].f_n[j]->max == 3)
         continue;
       for (k = 0; k < 4; k++)
         {
         nstore[0] += finalmesh->node[finalmesh->element[i].f_n[j]->list[k]].vert[0];
         nstore[1] += finalmesh->node[finalmesh->element[i].f_n[j]->list[k]].vert[1];
         nstore[2] += finalmesh->node[finalmesh->element[i].f_n[j]->list[k]].vert[2];
         } 

       fcent = Point(nstore[0]/4.0,nstore[1]/4.0,nstore[2]/4.0);

       normv = Vector(cent,fcent);
       normv.normalize();
       normv *= (fl*0.25);

       fcent += Point(normv[0],normv[1],normv[2]);

       //this is what we need to search through octree-sorted element list for!
       //set up pt, reset slop
       for (k = 0; k < 3; k++)
         {
         ptt[k] = fcent[k]; 
         stol1[k] = -1.0e+20;
         }
    
       //set up slop factor
       for (k = 0; k < 3; k++)
         stol1[k] = MAX(stol1[k],0.5*fl);
    
       //now, receive list of facets in the neighborhood to cut
       Octree_root1->retrieve_list(ptt,stol1,l,elist1);
    
       flag = 0; //jump out once found, check that is found
       for (k = 0; k < elist1->max && !flag; k++)
         {
         //set pntr pulled out of elist to what it actually is
         pntr1 = (POLY_ELEMENT*)elist1->list[k];

         //look to see if inside
         for (m = 0; m < pntr1->nf && !flag; m++)
           {
           if (pntr1->f_n[m]->max == 3)
             flag = line_facet_intersect(cent, fcent, finalmesh->node[pntr1->f_n[m]->list[0]].vert, finalmesh->node[pntr1->f_n[m]->list[1]].vert, finalmesh->node[pntr1->f_n[m]->list[2]].vert, pv, 0.1*fl);
           else
             {
             flag = line_facet_intersect(cent, fcent, finalmesh->node[pntr1->f_n[m]->list[2]].vert, finalmesh->node[pntr1->f_n[m]->list[1]].vert, finalmesh->node[pntr1->f_n[m]->list[0]].vert, pv, 0.1*fl);
             if (!flag)
               flag = line_facet_intersect(cent, fcent, finalmesh->node[pntr1->f_n[m]->list[3]].vert, finalmesh->node[pntr1->f_n[m]->list[2]].vert, finalmesh->node[pntr1->f_n[m]->list[0]].vert, pv, 0.1*fl);
             }
           }
         }

       if (!flag)
         {
         //set up diagonalized face and put in boundary list
         //clear tnode
         tnode->max = 0;
         for (k = 0; k < 4; k++)
           tnode->Add_To_List(finalmesh->element[i].f_n[j]->list[k]);

         tnt = 0;
         //call trimesh..will create elements and inc nelems
         finalmesh->set_up_tri_elem(tnt,tnode,map);

         //to keep track of trianlgles
         tempnt+= tnt;
         }
       }
     } 
    
   //finally, free octree mem 
   dummy1 = 0; 
   if (elist1 != 0)
     delete elist1;
   if (Octree_root1 != 0)
     delete Octree_root1;*/
   keepfacet.destruct();
   for (i = 0; i < geom->n_gnodes; i++)
     delete stack[i];
   delete stack;
   time(&tm);
   t_char = ctime(&tm);
   fprintf(out_f,"\nFinished viscous extrusion creation at %s",t_char);
   }
  
  /*//DEBUG
  for (i = 0; i < finalmesh->nn; i++)
    {
    for (j = i+1; j < finalmesh->nn; j++)
      {
      if (finalmesh->node[i].vert[0] == finalmesh->node[j].vert[0] && finalmesh->node[i].vert[1] == finalmesh->node[j].vert[1] && finalmesh->node[i].vert[2] == finalmesh->node[j].vert[2])
        {
        fprintf(out_f, "\nNode %d is the same as node %d\n",i,j);
        finalmesh->node[i].print(out_f);
        finalmesh->node[j].print(out_f);
        }
      }
    }*/
    
  //we will be setting up faces of internal voxels only based on neighbors, and we will need to add faces to voxels that are triangulated on the surface
  //this saves us being concerned about neighbors having same triangulation and needing to keep a list of tri faces per side
  //once all node numbers established, loop through and set up faces and surface mesh
  //we have all nodes defined, including mid-face and mid-edge, but NOT centroid of finest level
  //set up neighbor vars
  Point *pt = new Point[6];
  double begx, begy, begz, endx, endy, endz, halfx, halfy, halfz; 
  Point temp, pt0, pt1, pt2, pt3, pt4, pt5, pt6, pt7;
  Point *pnt = new Point[6];
  Vector *vec = new Vector[6];
  Vector pert;
  int nd[4], cnd[8], efnd[4][4], dfnd[4][6];
  tol = 0.25*smn;
  int fc = 0, q = 0, ind0 = 0;
  List dn;
  dn.construct();
  List **nbr = new List*[6];
  for (n = 0; n < 6; n++)
    nbr[n] = new List();

  flag = 0;
  //be sure tnt zeroed
  tnt = 0;
  
  if (!viscous)
    {
    tempnt = 0;
    //set up map for bd nodes
    map = new int[finalmesh->nn+geom->n_gnodes];
    //be sure init to -4
    for (k = 0; k < (finalmesh->nn+geom->n_gnodes); k++)
      map[k] = -4;

    //add in geom facets to tri
    ngt = 0;
  
    fprintf(out_f,"\nAllocating for geom boundary elements.\n");
    fflush(out_f);
  
    if ((finalmesh->nelems+geom->n_gfacets) > finalmesh->element_dim)
      my_mem(finalmesh->element_dim, finalmesh->nelems, &(finalmesh->element), (finalmesh->nelems+geom->n_gfacets), ELEMENT_CHUNK);
  
    //since all elements are one face, set up nf and f_n now
    for (k = finalmesh->nelems; k < (finalmesh->nelems+geom->n_gfacets); k++)
      {
      finalmesh->element[k].Initialize();
      finalmesh->element[k].nf = 1;
      finalmesh->element[k].f_n = new List*[1];
      finalmesh->element[k].f_n[0] = new List();
      }
    
    fprintf(out_f,"\nDetermining geometry tolerance and storing facets as boundary objects.\n");
    fflush(out_f);
      
    //reset tol for MIN op
    gtol = 1.0e+20;
  
    //loop through geom facets, noting all unique nodes at this point
    //set up geom as elements, then tris as elements
    for (i=1; i <= geom->ngb; i++)
      {
      for (j=geom->n_begin[i]; j < geom->n_end[i]; j++)
        {
        n0 = geom->g_facet[j].nodes[0];
        n1 = geom->g_facet[j].nodes[1];
        n2 = geom->g_facet[j].nodes[2];
      
        gtol = MIN(gtol,distance(geom->g_vert[n0], geom->g_vert[n1]));
        gtol = MIN(gtol,distance(geom->g_vert[n1], geom->g_vert[n2]));
        gtol = MIN(gtol,distance(geom->g_vert[n2], geom->g_vert[n0]));
      
        ngt++;
        //now, need to create tris
        if ((finalmesh->nelems+1) > finalmesh->element_dim)
          {
          my_mem(finalmesh->element_dim, finalmesh->nelems, &(finalmesh->element), (1+finalmesh->nelems), ELEMENT_CHUNK);
      
          //init elems
          for (n = finalmesh->nelems; n < (finalmesh->nelems+1); n++)
            {
            finalmesh->element[n].Initialize();
            finalmesh->element[n].nf = 1;
            finalmesh->element[n].f_n = new List*[1];
            finalmesh->element[n].f_n[0] = new List();
            }
          }  

        //take the resulting triangles and set up tri
        finalmesh->element[finalmesh->nelems].f_n[0]->Add_To_List(n0+finalmesh->nn);
        finalmesh->element[finalmesh->nelems].f_n[0]->Add_To_List(n1+finalmesh->nn);
        finalmesh->element[finalmesh->nelems].f_n[0]->Add_To_List(n2+finalmesh->nn);

        map[n0+finalmesh->nn] = -1;
        map[n1+finalmesh->nn] = -1;
        map[n2+finalmesh->nn] = -1;

      
        finalmesh->boundary[i-1].elist->Add_To_List(finalmesh->nelems); 
        //inc nelems
        finalmesh->nelems++;
        }
      }

    fprintf(out_f,"\n#Elem = %d\n",finalmesh->nelems);
    fflush(out_f); 

    //finalize gtol, print out so can be sent into restart
    gtol *= 0.25;
  
    fprintf(out_f,"\nGeometry tolerance is %16.10e.  Please enter after restart.\n",gtol);
    fflush(out_f);
  
    //set up tol as 0.25 smn, print out so can be sent into restart
    tol = 0.25*smn;
  
    fprintf(out_f,"\nInternal boundary tolerance is %16.10e.  Please enter after restart.\n",tol);
    fflush(out_f);
  
    //add in geom nodes so map is right
    for (i = 0; i < geom->n_gnodes; i++)
      {
      if ((finalmesh->nn + 1) > finalmesh->node_dim)
        my_mem(finalmesh->node_dim, finalmesh->nn, &(finalmesh->node), (finalmesh->nn + 1), NODE_CHUNK);
      
      finalmesh->node[finalmesh->nn].vert = geom->g_vert[i];
    
      //reset nn
      finalmesh->nn++;
      }
    
    if (ngt != geom->n_gfacets)
      {
      fprintf(out_f,"\nNumber of geometry facets doesn't match number to be written to UCD file. Exiting....\n");
      fflush(out_f);
      exit(0);
      }
    }
    
  time(&tm);
  t_char = ctime(&tm);
  fprintf(out_f,"\nTriangulating faces, allocating for internal elements, allocating for polymesh elements at %s.\n",t_char);
  fflush(out_f);  

  int ond[4];
    
  //start at j = 0 for elem numbering
  j = 0;
  for (k = 0; k < otmesh->nvoxels; k++)
    {
    if (otmesh->voxels[k].children.max > 0 || otmesh->voxels[k].cut != 1)
      continue; //kick out if children or if outside or if cut (no need for two directional check, since max 2 kids w/ 2 kids w/ orient known)
        
    //get crnerpts
    p0 = otmesh->voxels[k].cornerpts[0];
    p1 = otmesh->voxels[k].cornerpts[1];
  
    //get centroid
    cent = p0 + p1;
    cent /= 2.0;
  
    //set up values of other pts from cornerpts
    begx = p0[0];
    begy = p0[1];
    begz = p0[2];
    endx = p1[0];
    endy = p1[1];
    endz = p1[2];
    halfx = (p0[0]+p1[0])/2.0;
    halfy = (p0[1]+p1[1])/2.0; 
    halfz = (p0[2]+p1[2])/2.0;
          
    //set up other pts of this voxel
    pt[0] = Point(halfx,halfy,begz);
    pt[1] = Point(halfx,begy,halfz);
    pt[2] = Point(endx,halfy,halfz);
    pt[3] = Point(halfx,endy,halfz);
    pt[4] = Point(begx,halfy,halfz); 
    pt[5] = Point(halfx,halfy,endz);
      
    //set up cornerpts for use in face nabors
    pt0 = p0;
    pt1 = Point(endx,begy,begz);
    pt2 = Point(endx,endy,begz);
    pt3 = Point(begx,endy,begz);
    pt4 = Point(begx,begy,endz);
    pt5 = Point(endx,begy,endz);
    pt6 = p1;
    pt7 = Point(begx,endy,endz);

    //reset nbr
    for (n = 0; n < 6; n++)
      nbr[n]->max = 0;
      
    //be sure to reset elem face counter since constructor not called
    finalmesh->element[j].nf = 0;
    
    //now, we need to process and look at nabors...if cut surface, triangulate only cut part, if not, be sure we don't need to add quad faces
    //note...neighbor should ALWAYS exist!
    for (i = 0; i < 6; i++)
      {
      vec[i] = Vector(cent,pt[i]); //get direction
      vec[i].normalize(); //normalize
      vec[i] *= tol; //get tol as length
      pnt[i] = pt[i] + Point(vec[i][0],vec[i][1],vec[i][2]);
      
      //reset fc
      fc = 0;
      
      //now, branch out pnt four directions
      if (i == 0)
        {
        //set mid-edge nodes
        nd[0] = 16;
        nd[1] = 17;
        nd[2] = 18;
        nd[3] = 19;
        ond[0] = 8;
        ond[1] = 9;
        ond[2] = 10;
        ond[3] = 11;
        //set perts
        for (n = 0; n < 4; n++)
          {
          switch (n) 
            {
            case 0:
              pert = Vector(pt[i],pt0);
            break;
            case 1:
              pert = Vector(pt[i],pt1);
            break;
            case 2:
              pert = Vector(pt[i],pt3);
            break;
            case 3:
              pert = Vector(pt[i],pt2);
            break;
            default:
            break;
            }
          pert.normalize();
          pert *= tol;
          temp = pnt[i] + Point(pert[0],pert[1],pert[2]);
          nbr[i]->Add_To_List(otmesh->nabor_search(temp,0,smn,smx));
          }
        }
      else if (i == 1)
        {
        //set mid-edge nodes
        nd[0] = 10;
        nd[1] = 14;
        nd[2] = 15;
        nd[3] = 18;
        ond[0] = 8;
        ond[1] = 12;
        ond[2] = 13;
        ond[3] = 16;
        //set perts
        for (n = 0; n < 4; n++)
          {
          switch (n) 
            {
            case 0:
              pert = Vector(pt[i],pt0);
            break;
            case 1:
              pert = Vector(pt[i],pt1);
            break;
            case 2:
              pert = Vector(pt[i],pt4);
            break;
            case 3:
              pert = Vector(pt[i],pt5);
            break;
            default:
            break;
            }
          pert.normalize();
          pert *= tol;
          temp = pnt[i] + Point(pert[0],pert[1],pert[2]);
          nbr[i]->Add_To_List(otmesh->nabor_search(temp,0,smn,smx));
          }
        }
      else if (i == 2)
        {
        //set mid-edge nodes
        nd[0] = 11;
        nd[1] = 12;
        nd[2] = 15;
        nd[3] = 19;
        ond[0] = 9;
        ond[1] = 13;
        ond[2] = 14;
        ond[3] = 17;
        //set perts
        for (n = 0; n < 4; n++)
          {
          switch (n) 
            {
            case 0:
              pert = Vector(pt[i],pt1);
            break;
            case 1:
              pert = Vector(pt[i],pt2);
            break;
            case 2:
              pert = Vector(pt[i],pt5);
            break;
            case 3:
              pert = Vector(pt[i],pt6);
            break;
            default:
            break;
            }
          pert.normalize();
          pert *= tol;
          temp = pnt[i] + Point(pert[0],pert[1],pert[2]);
          nbr[i]->Add_To_List(otmesh->nabor_search(temp,0,smn,smx));
          }
        }
      else if (i == 3)
        {
        //set mid-edge nodes
        nd[0] = 8;
        nd[1] = 12;
        nd[2] = 13;
        nd[3] = 16;
        ond[0] = 10;
        ond[1] = 14;
        ond[2] = 15;
        ond[3] = 18;
        //set perts
        for (n = 0; n < 4; n++)
          {
          switch (n) 
            {
            case 0:
              pert = Vector(pt[i],pt3);
            break;
            case 1:
              pert = Vector(pt[i],pt2);
            break;
            case 2:
              pert = Vector(pt[i],pt7);
            break;
            case 3:
              pert = Vector(pt[i],pt6);
            break;
            default:
            break;
            }
          pert.normalize();
          pert *= tol;
          temp = pnt[i] + Point(pert[0],pert[1],pert[2]);
          nbr[i]->Add_To_List(otmesh->nabor_search(temp,0,smn,smx));
          }
        }
      else if (i == 4)
        {
        //set mid-edge nodes
        nd[0] = 9;
        nd[1] = 13;
        nd[2] = 14;
        nd[3] = 17;
        ond[0] = 11;
        ond[1] = 12;
        ond[2] = 15;
        ond[3] = 19;
        //set perts
        for (n = 0; n < 4; n++)
          {
          switch (n) 
            {
            case 0:
              pert = Vector(pt[i],pt0);
            break;
            case 1:
              pert = Vector(pt[i],pt3);
            break;
            case 2:
              pert = Vector(pt[i],pt4);
            break;
            case 3:
              pert = Vector(pt[i],pt7);
            break;
            default:
            break;
            }
          pert.normalize();
          pert *= tol;
          temp = pnt[i] + Point(pert[0],pert[1],pert[2]);
          nbr[i]->Add_To_List(otmesh->nabor_search(temp,0,smn,smx));
          }
        }
      else if (i == 5)
        {
        //set mid-edge nodes
        nd[0] = 8;
        nd[1] = 9;
        nd[2] = 10;
        nd[3] = 11;
        ond[0] = 16;
        ond[1] = 17;
        ond[2] = 18;
        ond[3] = 19;
        //set perts
        for (n = 0; n < 4; n++)
          {
          switch (n) 
            {
            case 0:
              pert = Vector(pt[i],pt4);
            break;
            case 1:
              pert = Vector(pt[i],pt5);
            break;
            case 2:
              pert = Vector(pt[i],pt7);
            break;
            case 3:
              pert = Vector(pt[i],pt6);
            break;
            default:
            break;
            }
          pert.normalize();
          pert *= tol;
          temp = pnt[i] + Point(pert[0],pert[1],pert[2]);
          nbr[i]->Add_To_List(otmesh->nabor_search(temp,0,smn,smx));
          }
        }
      
      int oversetflag = 0;
      if (overset)
        {
        if (nbr[i]->Times_In_List(nbr[i]->list[0]) == 4 && nbr[i]->list[0] == -1)
          {
          for (m = 0; m < 4; m++)
            if (otmesh->voxels[k].nodes[ond[m]] >= 0)
              fc++;
          finalmesh->element[j].nf += (fc+2);
          oversetflag++;
          }
        }

      //set up list to hold done neighbors
      dn.Redimension(0);
      flag = 0;
      //now, increment face counts based on cut as well as num neighbors
      for (n = 0; n < 4 && !flag && !oversetflag; n++)
        {
        if (nbr[i]->Times_In_List(nbr[i]->list[n]) == 4)
          {
          if (otmesh->voxels[nbr[i]->list[n]].cut == 2)
            {
            for (m = 0; m < 4; m++)
              if (otmesh->voxels[nbr[i]->list[n]].nodes[nd[m]] >= 0)
                fc++;
            finalmesh->element[j].nf += (fc+2);    
            }
          else
            finalmesh->element[j].nf++; //only have one nabor and only need one face
          //flag and hop out of loop  
          flag++;
          }
        //nothing will ever be in list 3 times...
        else if (nbr[i]->Times_In_List(nbr[i]->list[n]) == 2 && !dn.Is_In_List(nbr[i]->list[n]))
          {
          if (otmesh->voxels[nbr[i]->list[n]].cut == 2)
            {
            for (m = 0; m < 4; m++)
              if (otmesh->voxels[nbr[i]->list[n]].nodes[nd[m]] >= 0)
                fc++;
            finalmesh->element[j].nf += (fc+2);    
            }
          else
            finalmesh->element[j].nf++; //not triangulated, so just one face
          //add to done list
          dn.Add_To_List(nbr[i]->list[n]);
          }
        else if (nbr[i]->Times_In_List(nbr[i]->list[n]) == 1)
          {
          if (otmesh->voxels[nbr[i]->list[n]].cut == 2)
            {
            finalmesh->element[j].nf += 2; //this is guaranteed to be a quarter face, and by qual const (4-1 caddy corner), it cannot have more than four nodes to triangulate   
            }
          else
            finalmesh->element[j].nf++; //not triangulated, so just one face
          }
        }
      }
        
    //now, set up proper lists
    finalmesh->element[j].f_n = new List*[finalmesh->element[j].nf];
    for (i = 0; i < finalmesh->element[j].nf; i++)
      finalmesh->element[j].f_n[i] = new List();
    
    //reset fc for face counter
    fc = 0;
    
    //now, create cut surface and element faces side by side using nbr already created        
    for (i = 0; i < 6; i++)
      {
      switch (i) 
        {
        case 0:
          cnd[0] = 0;
          cnd[1] = 11;
          cnd[2] = 3;
          cnd[3] = 10;
          cnd[4] = 2;
          cnd[5] = 9;
          cnd[6] = 1;
          cnd[7] = 8;
          efnd[0][0] = 11;
          efnd[0][1] = 20;
          efnd[0][2] = 8;
          efnd[0][3] = 0;
          efnd[1][0] = 20;
          efnd[1][1] = 9;
          efnd[1][2] = 1;
          efnd[1][3] = 8;
          efnd[2][0] = 3;
          efnd[2][1] = 10;
          efnd[2][2] = 20;
          efnd[2][3] = 11;
          efnd[3][0] = 10;
          efnd[3][1] = 2;
          efnd[3][2] = 9;
          efnd[3][3] = 20;
          dfnd[0][0] = 3;
          dfnd[0][1] = 10;
          dfnd[0][2] = 20;
          dfnd[0][3] = 8;
          dfnd[0][4] = 0;
          dfnd[0][5] = 11;
          dfnd[1][0] = 10;
          dfnd[1][1] = 2;
          dfnd[1][2] = 9;
          dfnd[1][3] = 1;
          dfnd[1][4] = 8;
          dfnd[1][5] = 20;
          dfnd[2][0] = 20;
          dfnd[2][1] = 9;
          dfnd[2][2] = 1;
          dfnd[2][3] = 8;
          dfnd[2][4] = 0;
          dfnd[2][5] = 11;
          dfnd[3][0] = 3;
          dfnd[3][1] = 10;
          dfnd[3][2] = 2;
          dfnd[3][3] = 9;
          dfnd[3][4] = 20;
          dfnd[3][5] = 11;
        break;
        case 1:
          cnd[0] = 0;
          cnd[1] = 8;
          cnd[2] = 1;
          cnd[3] = 13;
          cnd[4] = 5;
          cnd[5] = 16;
          cnd[6] = 4;
          cnd[7] = 12;
          efnd[0][0] = 0;
          efnd[0][1] = 8;
          efnd[0][2] = 21;
          efnd[0][3] = 12;
          efnd[1][0] = 8;
          efnd[1][1] = 1;
          efnd[1][2] = 13;
          efnd[1][3] = 21;
          efnd[2][0] = 12;
          efnd[2][1] = 21;
          efnd[2][2] = 16;
          efnd[2][3] = 4;
          efnd[3][0] = 21;
          efnd[3][1] = 13;
          efnd[3][2] = 5;
          efnd[3][3] = 16;
          dfnd[0][0] = 0;
          dfnd[0][1] = 8;
          dfnd[0][2] = 21;
          dfnd[0][3] = 16;
          dfnd[0][4] = 4;
          dfnd[0][5] = 12;
          dfnd[1][0] = 8;
          dfnd[1][1] = 1;
          dfnd[1][2] = 13;
          dfnd[1][3] = 5;
          dfnd[1][4] = 16;
          dfnd[1][5] = 21;
          dfnd[2][0] = 0;
          dfnd[2][1] = 8;
          dfnd[2][2] = 1;
          dfnd[2][3] = 13;
          dfnd[2][4] = 21;
          dfnd[2][5] = 12;
          dfnd[3][0] = 12;
          dfnd[3][1] = 21;
          dfnd[3][2] = 13;
          dfnd[3][3] = 5;
          dfnd[3][4] = 16;
          dfnd[3][5] = 4;
        break;
        case 2:
          cnd[0] = 1;
          cnd[1] = 9;
          cnd[2] = 2;
          cnd[3] = 14;
          cnd[4] = 6;
          cnd[5] = 17;
          cnd[6] = 5;
          cnd[7] = 13;
          efnd[0][0] = 1;
          efnd[0][1] = 9;
          efnd[0][2] = 22;
          efnd[0][3] = 13;
          efnd[1][0] = 9;
          efnd[1][1] = 2;
          efnd[1][2] = 14;
          efnd[1][3] = 22;
          efnd[2][0] = 13;
          efnd[2][1] = 22;
          efnd[2][2] = 17;
          efnd[2][3] = 5;
          efnd[3][0] = 22;
          efnd[3][1] = 14;
          efnd[3][2] = 6;
          efnd[3][3] = 17;
          dfnd[0][0] = 1;
          dfnd[0][1] = 9;
          dfnd[0][2] = 22;
          dfnd[0][3] = 17;
          dfnd[0][4] = 5;
          dfnd[0][5] = 13;
          dfnd[1][0] = 9;
          dfnd[1][1] = 2;
          dfnd[1][2] = 14;
          dfnd[1][3] = 6;
          dfnd[1][4] = 17;
          dfnd[1][5] = 22;
          dfnd[2][0] = 1;
          dfnd[2][1] = 9;
          dfnd[2][2] = 2;
          dfnd[2][3] = 14;
          dfnd[2][4] = 22;
          dfnd[2][5] = 13;
          dfnd[3][0] = 13;
          dfnd[3][1] = 22;
          dfnd[3][2] = 14;
          dfnd[3][3] = 6;
          dfnd[3][4] = 17;
          dfnd[3][5] = 5;
        break;
        case 3:
          cnd[0] = 2;
          cnd[1] = 10;
          cnd[2] = 3;
          cnd[3] = 15;
          cnd[4] = 7;
          cnd[5] = 18;
          cnd[6] = 6;
          cnd[7] = 14;
          efnd[0][0] = 10;
          efnd[0][1] = 3;
          efnd[0][2] = 15;
          efnd[0][3] = 23;
          efnd[1][0] = 2;
          efnd[1][1] = 10;
          efnd[1][2] = 23;
          efnd[1][3] = 14;
          efnd[2][0] = 23;
          efnd[2][1] = 15;
          efnd[2][2] = 7;
          efnd[2][3] = 18;
          efnd[3][0] = 14;
          efnd[3][1] = 23;
          efnd[3][2] = 18;
          efnd[3][3] = 6;
          dfnd[0][0] = 10;
          dfnd[0][1] = 3;
          dfnd[0][2] = 15;
          dfnd[0][3] = 7;
          dfnd[0][4] = 18;
          dfnd[0][5] = 23;
          dfnd[1][0] = 2;
          dfnd[1][1] = 10;
          dfnd[1][2] = 23;
          dfnd[1][3] = 18;
          dfnd[1][4] = 6;
          dfnd[1][5] = 14;
          dfnd[2][0] = 2;
          dfnd[2][1] = 10;
          dfnd[2][2] = 3;
          dfnd[2][3] = 15;
          dfnd[2][4] = 23;
          dfnd[2][5] = 14;
          dfnd[3][0] = 14;
          dfnd[3][1] = 23;
          dfnd[3][2] = 15;
          dfnd[3][3] = 7;
          dfnd[3][4] = 18;
          dfnd[3][5] = 6;
        break;
        case 4:
          cnd[0] = 3;
          cnd[1] = 11;
          cnd[2] = 0;
          cnd[3] = 12;
          cnd[4] = 4;
          cnd[5] = 19;
          cnd[6] = 7;
          cnd[7] = 15;
          efnd[0][0] = 11;
          efnd[0][1] = 0;
          efnd[0][2] = 12;
          efnd[0][3] = 24;
          efnd[1][0] = 3;
          efnd[1][1] = 11;
          efnd[1][2] = 24;
          efnd[1][3] = 15;
          efnd[2][0] = 24;
          efnd[2][1] = 12;
          efnd[2][2] = 4;
          efnd[2][3] = 19;
          efnd[3][0] = 15;
          efnd[3][1] = 24;
          efnd[3][2] = 19;
          efnd[3][3] = 7;
          dfnd[0][0] = 11;
          dfnd[0][1] = 0;
          dfnd[0][2] = 12;
          dfnd[0][3] = 4;
          dfnd[0][4] = 19;
          dfnd[0][5] = 24;
          dfnd[1][0] = 3;
          dfnd[1][1] = 11;
          dfnd[1][2] = 24;
          dfnd[1][3] = 19;
          dfnd[1][4] = 7;
          dfnd[1][5] = 15;
          dfnd[2][0] = 3;
          dfnd[2][1] = 11;
          dfnd[2][2] = 0;
          dfnd[2][3] = 12;
          dfnd[2][4] = 24;
          dfnd[2][5] = 15;
          dfnd[3][0] = 15;
          dfnd[3][1] = 24;
          dfnd[3][2] = 12;
          dfnd[3][3] = 4;
          dfnd[3][4] = 19;
          dfnd[3][5] = 7;
        break;
        case 5:
          cnd[0] = 6;
          cnd[1] = 18;
          cnd[2] = 7;
          cnd[3] = 19;
          cnd[4] = 4;
          cnd[5] = 16;
          cnd[6] = 5;
          cnd[7] = 17;
          efnd[0][0] = 25;
          efnd[0][1] = 19;
          efnd[0][2] = 4;
          efnd[0][3] = 16;
          efnd[1][0] = 17;
          efnd[1][1] = 25;
          efnd[1][2] = 16;
          efnd[1][3] = 5;
          efnd[2][0] = 18;
          efnd[2][1] = 7;
          efnd[2][2] = 19;
          efnd[2][3] = 25;
          efnd[3][0] = 6;
          efnd[3][1] = 18;
          efnd[3][2] = 25;
          efnd[3][3] = 17;
          dfnd[0][0] = 18;
          dfnd[0][1] = 7;
          dfnd[0][2] = 19;
          dfnd[0][3] = 4;
          dfnd[0][4] = 16;
          dfnd[0][5] = 25;
          dfnd[1][0] = 6;
          dfnd[1][1] = 18;
          dfnd[1][2] = 25;
          dfnd[1][3] = 16;
          dfnd[1][4] = 5;
          dfnd[1][5] = 17;
          dfnd[2][0] = 17;
          dfnd[2][1] = 25;
          dfnd[2][2] = 19;
          dfnd[2][3] = 4;
          dfnd[2][4] = 16;
          dfnd[2][5] = 5;
          dfnd[3][0] = 6;
          dfnd[3][1] = 18;
          dfnd[3][2] = 7;
          dfnd[3][3] = 19;
          dfnd[3][4] = 25;
          dfnd[3][5] = 17;
        break;
        default:
        break;
        }
      
      int oversetflag = 0;
      if (overset)
        {
        //clear tnode
        tnode->max = 0;
        if (nbr[i]->Times_In_List(nbr[i]->list[0]) == 4 && nbr[i]->list[0] == -1)
          {
          for (m = 0; m < 8; m++)
            if (otmesh->voxels[k].nodes[cnd[m]] >= 0)
              tnode->Add_To_List(otmesh->voxels[k].nodes[cnd[m]]);
            
          tnt = 0;
          //call trimesh    

          //tnode->print(out_f);

          finalmesh->set_up_tri(tnt,tnode,map);
            
          //now, make these tri into faces
          for (m = 0; m < tnt; m++)
            {
            for (q = 0; q < finalmesh->element[finalmesh->nelems + m].f_n[0]->max; q++)
              finalmesh->element[j].f_n[fc]->Add_To_List(finalmesh->element[finalmesh->nelems + m].f_n[0]->list[q]);
            fc++; //added face
            }
          //now, inc nelems
          finalmesh->nelems+= tnt;
          //to keep track of trianlgles
          tempnt+= tnt;
          oversetflag++;
          }
        }

      //set up list to hold done neighbors
      dn.Redimension(0);
      flag = 0;
      //now, increment face counts based on cut as well as num neighbors
      for (n = 0; n < 4 && !flag && !oversetflag; n++)
        {
        //clear tnode
        tnode->max = 0;
        //look at sole nabor
        if (nbr[i]->Times_In_List(nbr[i]->list[n]) == 4)
          {
          if (otmesh->voxels[nbr[i]->list[n]].cut == 2)
            {
            for (m = 0; m < 8; m++)
              if (otmesh->voxels[k].nodes[cnd[m]] >= 0)
                tnode->Add_To_List(otmesh->voxels[k].nodes[cnd[m]]);
            
            tnt = 0;
            //call trimesh    
            finalmesh->set_up_tri(tnt,tnode,map);
            
            //now, make these tri into faces
            for (m = 0; m < tnt; m++)
              {
              for (q = 0; q < finalmesh->element[finalmesh->nelems + m].f_n[0]->max; q++)
                finalmesh->element[j].f_n[fc]->Add_To_List(finalmesh->element[finalmesh->nelems + m].f_n[0]->list[q]);
              fc++; //added face
              }
            //now, inc nelems
            finalmesh->nelems+= tnt;
            //to keep track of trianlgles
            tempnt+= tnt;
            }
          else
            {
            for (m = 0; m < 8; m++)
              if (otmesh->voxels[k].nodes[cnd[m]] >= 0)
                finalmesh->element[j].f_n[fc]->Add_To_List(otmesh->voxels[k].nodes[cnd[m]]);
              
            fc++; //added face
            }
          //flag and hop out of loop  
          flag++;
          }
        //nothing will ever be in list 3 times...
        else if (nbr[i]->Times_In_List(nbr[i]->list[n]) == 2 && !dn.Is_In_List(nbr[i]->list[n]))
          {
          if (otmesh->voxels[nbr[i]->list[n]].cut == 2)
            {
            //need to determine which two nabors this is
            for (q = 0; q < 4; q++)
              {
              if (q == n)
                continue;
              if (nbr[i]->list[n] == nbr[i]->list[q])
                ind0 = q; 
              }
            if ((ind0 == 0 && n == 2) || (ind0 == 2 && n == 0))
              q = 0;
            if ((ind0 == 1 && n == 3) || (ind0 == 3 && n == 1))
              q = 1;
            if ((ind0 == 0 && n == 1) || (ind0 == 1 && n == 0))
              q = 2;
            if ((ind0 == 2 && n == 3) || (ind0 == 3 && n == 2))
              q = 3;
            //now set up faces
            for (m = 0; m < 6; m++)
              if (otmesh->voxels[k].nodes[dfnd[q][m]] >= 0)
                tnode->Add_To_List(otmesh->voxels[k].nodes[dfnd[q][m]]);
              
            tnt = 0;
            //call trimesh    
            finalmesh->set_up_tri(tnt,tnode,map);
            
            //now, make these tri into faces
            for (m = 0; m < tnt; m++)
              {
              for (q = 0; q < finalmesh->element[finalmesh->nelems + m].f_n[0]->max; q++)
                finalmesh->element[j].f_n[fc]->Add_To_List(finalmesh->element[finalmesh->nelems + m].f_n[0]->list[q]);
              fc++; //added face
              }
            //now, inc nelems
            finalmesh->nelems+= tnt;
            //to keep track of trianlgles
            tempnt+= tnt;
            }
          else
            {
            //need to determine which two nabors this is
            for (q = 0; q < 4; q++)
              {
              if (q == n)
                continue;
              if (nbr[i]->list[n] == nbr[i]->list[q])
                ind0 = q; 
              }
            if ((ind0 == 0 && n == 2) || (ind0 == 2 && n == 0))
              q = 0;
            if ((ind0 == 1 && n == 3) || (ind0 == 3 && n == 1))
              q = 1;
            if ((ind0 == 0 && n == 1) || (ind0 == 1 && n == 0))
              q = 2;
            if ((ind0 == 2 && n == 3) || (ind0 == 3 && n == 2))
              q = 3;
            //now set up faces
            for (m = 0; m < 6; m++)
              if (otmesh->voxels[k].nodes[dfnd[q][m]] >= 0)
                finalmesh->element[j].f_n[fc]->Add_To_List(otmesh->voxels[k].nodes[dfnd[q][m]]);
            fc++;
            }
          //add to done list
          dn.Add_To_List(nbr[i]->list[n]);
          }
        else if (nbr[i]->Times_In_List(nbr[i]->list[n]) == 1)
          {
          if (otmesh->voxels[nbr[i]->list[n]].cut == 2)
            {
            for (m = 0; m < 4; m++)
              {
              if (otmesh->voxels[k].nodes[efnd[n][m]] < 0)
                {
                fprintf(stderr, "\nYou have an undefined node in voxel %d at position %d, m= %d, n= %d, i = %d.  Exiting....\n",k,efnd[n][m],m,n,i);
                otmesh->voxels[k].print(stderr);
                nbr[i]->print(stderr);
                otmesh->voxels[nbr[i]->list[0]].print(stderr);
                otmesh->voxels[nbr[i]->list[1]].print(stderr);
                otmesh->voxels[nbr[i]->list[2]].print(stderr);
                otmesh->voxels[nbr[i]->list[3]].print(stderr);
                fflush(stderr);
                exit(0);
                }
              tnode->Add_To_List(otmesh->voxels[k].nodes[efnd[n][m]]);
              }
            
            tnt = 0;
            //call trimesh    
            finalmesh->set_up_tri(tnt,tnode,map);
            
            //now, make these tri into faces
            for (m = 0; m < tnt; m++)
              {
              for (q = 0; q < finalmesh->element[finalmesh->nelems + m].f_n[0]->max; q++)
                finalmesh->element[j].f_n[fc]->Add_To_List(finalmesh->element[finalmesh->nelems + m].f_n[0]->list[q]);
              fc++; //added face
              }
            //now, inc nelems
            finalmesh->nelems+= tnt;
            //to keep track of trianlgles
            tempnt+= tnt;
            }
          else
            {
            for (m = 0; m < 4; m++)
              {
              if (otmesh->voxels[k].nodes[efnd[n][m]] < 0)
                {
                fprintf(stderr, "\nYou have an undefined node in voxel %d at position %d, m= %d, n= %d, i = %d.  Exiting....\n",k,efnd[n][m],m,n,i);
                fflush(stderr);
                otmesh->voxels[k].print(stderr);
                nbr[i]->print(stderr);
                otmesh->voxels[nbr[i]->list[0]].print(stderr);
                otmesh->voxels[nbr[i]->list[1]].print(stderr);
                otmesh->voxels[nbr[i]->list[2]].print(stderr);
                otmesh->voxels[nbr[i]->list[3]].print(stderr);
                fflush(stderr);
                exit(0);
                }
              finalmesh->element[j].f_n[fc]->Add_To_List(otmesh->voxels[k].nodes[efnd[n][m]]);
              }
            fc++;
            }
          }
        }
      }
    j++; //inc element
    }
    
  //delete nbr
  for (i = 0; i < 6; i++)
    delete nbr[i];
  delete [] nbr;
  delete tnode;
  dn.destruct();
  
  //delete nabor info
  delete [] vec;
  delete [] pt;
  delete [] pnt;
  
  //VCB: FOR NOW, DELETE OTMESH..no tree needed until we go dynamic...save space
  delete otmesh;

  //now, cut things off if overset
  if (overset)
    {
    //reset file name before volume mesh output
    sname[0] = '\0';
  
    sprintf(sname,"FASTAR_cartesian.cgns");

    //output vol mesh
    finalmesh->smooth_io(1,geom,sname);

    //reset file name before internal mesh output
    sname[0] = '\0';
  
    sprintf(sname,"FASTAR_geom.cgns");

    //output vol mesh
    ovst->smooth_io(1,geom,sname);
  
    fprintf(out_f,"\nFinished with output CGNS files.\n");
    fflush(out_f);

    delete finalmesh;
    delete ovst;

    return;
    }
  
  fprintf(out_f,"\nDone creating polymesh. Readying map for UCD/STL file.\n");
  fflush(out_f);
  
  //NOTE: Bruce's tet mesher will go here, and avoid the rest to //***, leaving restart in
  //Also, if we can pass him any voxel nodes that are in cut voxels but inside the computational domain as well, he can use them as seed nodes.
  //Just add to mesh obj as last nodes (will be kept anyhow) and give starting index
    
  snn = 0;
  gnn = 0;
  nni = 0;
    
  //determine number of nodes that will be added to surface file, and remap nodes
  //will mix up bd types, but not important for looping through elements and will come back out of order anyways
  for (i = 0; i < finalmesh->nn; i++)
    {
    if (map[i] == -2)
      {
      map[i] = nni++;
      snn++;
      }
    if (map[i] == -1)
      {
      map[i] = nni++;
      gnn++;
      }
    }
      
  fprintf(out_f,"\nNumber of nodes on internal stitching boundary = %d\n",snn);
  fflush(out_f);
  
  fprintf(out_f,"\nNumber of triangles on internal stitching boundary = %d\n",tempnt);
  fflush(out_f);
      
  fprintf(out_f,"\nNumber of nodes on geometry boundary = %d\n",gnn);
  fflush(out_f);
  
  fprintf(out_f,"\nNumber of triangles on geometry boundary = %d\n",ngt);
  fflush(out_f);
  
  //finally, output ucd file
  sname[0] = '\0';
  
  if (!ftype)
    sprintf(sname,"FASTAR.ucd");
  else
    sprintf(sname,"FASTAR.stl");

  if ((ucd_f=fopen(sname,"w")) == NULL)
    {
    fprintf(stderr,"\nCouldn't open surface mesh output file %s",sname);
    fflush(stderr);
    //MPI_Abort(MPI_COMM_WORLD,0);
    exit(0);
    }
   
  if (!ftype)
    {   
    //put in totals of elements and nodes
    if (!viscous)
      {
      fprintf(ucd_f, "%d %d %d %d %d\n",nni,tempnt+ngt,0,0,0);
      fflush(ucd_f);
      }
    else
      {
      fprintf(ucd_f, "%d %d %d %d %d\n",nni,tempnt+nonviscfacets,0,0,0);
      fflush(ucd_f);
      }

    //need to put in remapped finalmesh pts, then geom pts (already "remapped")
    for (i = 0; i < finalmesh->nn; i++)
      if (map[i] >= 0)
        fprintf(ucd_f, "%d %16.10e %16.10e %16.10e\n",map[i]+1,finalmesh->node[i].vert[0],finalmesh->node[i].vert[1],finalmesh->node[i].vert[2]);
    fflush(ucd_f);
    
    //need to put in surface tri and geom tri, remapping pts
    int tempin = 1;
    Vector areasum = Vector(0.0,0.0,0.0);
    int chkbd = geom->ngb+1;
    if (viscous)
      chkbd++;

    for (i = 0; i < chkbd; i++)
      {
      if (viscous)
        if (!nonviscbd.Is_In_List(i+1))
          continue;

      for (j = 0; j < finalmesh->boundary[i].elist->max; j++)
        {
        n0 = finalmesh->element[finalmesh->boundary[i].elist->list[j]].f_n[0]->list[0];
        n1 = finalmesh->element[finalmesh->boundary[i].elist->list[j]].f_n[0]->list[1];
        n2 = finalmesh->element[finalmesh->boundary[i].elist->list[j]].f_n[0]->list[2];
      
        //find vectors for area sum
        v0 = Vector(finalmesh->node[n0].vert,finalmesh->node[n1].vert);
        v1 = Vector(finalmesh->node[n0].vert,finalmesh->node[n2].vert);
        
        //do area summation for inner bd only!
        if (i == geom->ngb)
          {
          norm = (v0%v1)*0.5;
          areasum += norm;
          }
        
        fprintf(ucd_f, "%d %d tri %d %d %d\n",tempin,i,map[n0]+1,map[n1]+1,map[n2]+1);
        tempin++;
        }
      }
    fflush(ucd_f);

    fprintf(out_f, "Area summation (x, y, z) = %16.10e %16.10e %16.10e\n",areasum[0],areasum[1],areasum[2]);
    fflush(out_f);
    
    //now, delete map
    if (!gentet)
      delete [] map;
    }
  else if (ftype == -1)
    {
    //now, delete map
    if (!gentet)
      delete [] map;
    
    //STL format
    fprintf(ucd_f, "solid for_pointwise\n");
    fflush(ucd_f);
    
    Vector areasum = Vector(0.0,0.0,0.0);
    
    //put in triangles and normals
    int chkbd = geom->ngb+1;
    if (viscous)
      chkbd++;

    for (i = 0; i < chkbd; i++)
      {
      //puts viscous stitch in with geom
      if (viscous)
        if (!nonviscbd.Is_In_List(i+1) || i == geom->ngb)
          continue;

      for (j = 0; j < finalmesh->boundary[i].elist->max; j++)
        {
        n0 = finalmesh->element[finalmesh->boundary[i].elist->list[j]].f_n[0]->list[0];
        n1 = finalmesh->element[finalmesh->boundary[i].elist->list[j]].f_n[0]->list[1];
        n2 = finalmesh->element[finalmesh->boundary[i].elist->list[j]].f_n[0]->list[2];
      
        //find normal vector
        v0 = Vector(finalmesh->node[n0].vert,finalmesh->node[n1].vert);
        v1 = Vector(finalmesh->node[n0].vert,finalmesh->node[n2].vert);
        norm = v0 % v1;
        norm.normalize();
        
        fprintf(ucd_f, "facet normal %16.10e %16.10e %16.10e\n",norm[0],norm[1],norm[2]);
        
        fprintf(ucd_f, "  outer loop\n");
        fprintf(ucd_f, "    vertex %16.10e %16.10e %16.10e\n",finalmesh->node[n0].vert[0],finalmesh->node[n0].vert[1],finalmesh->node[n0].vert[2]);
        fprintf(ucd_f, "    vertex %16.10e %16.10e %16.10e\n",finalmesh->node[n1].vert[0],finalmesh->node[n1].vert[1],finalmesh->node[n1].vert[2]);
        fprintf(ucd_f, "    vertex %16.10e %16.10e %16.10e\n",finalmesh->node[n2].vert[0],finalmesh->node[n2].vert[1],finalmesh->node[n2].vert[2]);
        fprintf(ucd_f, "  endloop\n");
        fprintf(ucd_f, "endfacet\n");
        }
      }
      
    fprintf(ucd_f, "endsolid for_pointwise\n");
    fflush(ucd_f);

    //this is klugy...making two stl files for inner and outer...needless to say, with cube in cube, we'd have to separate stitching bds which would be bad!
    //close ucd file
    fclose(ucd_f);

    //reopen
    sname[0] = '\0';
    sprintf(sname,"FASTAR_inner.stl");

    if ((ucd_f=fopen(sname,"w")) == NULL)
      {
      fprintf(stderr,"\nCouldn't open surface mesh output file %s",sname);
      fflush(stderr);
      //MPI_Abort(MPI_COMM_WORLD,0);
      exit(0);
      }

    //STL format
    fprintf(ucd_f, "solid for_pointwise2\n");
    fflush(ucd_f);
    
    areasum = Vector(0.0,0.0,0.0);
    
    //put in triangles and normals
    for (i = geom->ngb; i < geom->ngb+1; i++)
      {
      for (j = 0; j < finalmesh->boundary[i].elist->max; j++)
        {
        n0 = finalmesh->element[finalmesh->boundary[i].elist->list[j]].f_n[0]->list[0];
        n1 = finalmesh->element[finalmesh->boundary[i].elist->list[j]].f_n[0]->list[1];
        n2 = finalmesh->element[finalmesh->boundary[i].elist->list[j]].f_n[0]->list[2];
      
        //find normal vector
        v0 = Vector(finalmesh->node[n0].vert,finalmesh->node[n1].vert);
        v1 = Vector(finalmesh->node[n0].vert,finalmesh->node[n2].vert);
        norm = v0 % v1;
        norm.normalize();
        
        fprintf(ucd_f, "facet normal %16.10e %16.10e %16.10e\n",norm[0],norm[1],norm[2]);
        
        //do area summation for inner bd only..reset norm in process
        norm = (v0%v1)*0.5;
        areasum += norm;
        
        fprintf(ucd_f, "  outer loop\n");
        fprintf(ucd_f, "    vertex %16.10e %16.10e %16.10e\n",finalmesh->node[n0].vert[0],finalmesh->node[n0].vert[1],finalmesh->node[n0].vert[2]);
        fprintf(ucd_f, "    vertex %16.10e %16.10e %16.10e\n",finalmesh->node[n1].vert[0],finalmesh->node[n1].vert[1],finalmesh->node[n1].vert[2]);
        fprintf(ucd_f, "    vertex %16.10e %16.10e %16.10e\n",finalmesh->node[n2].vert[0],finalmesh->node[n2].vert[1],finalmesh->node[n2].vert[2]);
        fprintf(ucd_f, "  endloop\n");
        fprintf(ucd_f, "endfacet\n");
        }
      }

    fprintf(ucd_f, "endsolid for_pointwise2\n");
    fflush(ucd_f);

    fprintf(out_f, "Area summation Two (x, y, z) = %16.10e %16.10e %16.10e\n",areasum[0],areasum[1],areasum[2]);
    fflush(out_f);
    }
  else if (ftype == 1)
    {
    //now, delete map
    if (!gentet)
      delete [] map;
    
    //STL format
    fprintf(ucd_f, "solid for_pointwise\n");
    fflush(ucd_f);
    
    Vector areasum = Vector(0.0,0.0,0.0);
    
    //put in triangles and normals
    int chkbd = geom->ngb+1;
    if (viscous)
      chkbd++;

    for (i = 0; i < chkbd; i++)
      {
      if (viscous)
        if (!nonviscbd.Is_In_List(i+1))
          continue;

      for (j = 0; j < finalmesh->boundary[i].elist->max; j++)
        {
        n0 = finalmesh->element[finalmesh->boundary[i].elist->list[j]].f_n[0]->list[0];
        n1 = finalmesh->element[finalmesh->boundary[i].elist->list[j]].f_n[0]->list[1];
        n2 = finalmesh->element[finalmesh->boundary[i].elist->list[j]].f_n[0]->list[2];
      
        //find normal vector
        v0 = Vector(finalmesh->node[n0].vert,finalmesh->node[n1].vert);
        v1 = Vector(finalmesh->node[n0].vert,finalmesh->node[n2].vert);
        norm = v0 % v1;
        norm.normalize();
        
        fprintf(ucd_f, "facet normal %16.10e %16.10e %16.10e\n",norm[0],norm[1],norm[2]);
        
        //do area summation for inner bd only..reset norm in process
        if (i == geom->ngb)
          {
          norm = (v0%v1)*0.5;
          areasum += norm;
          }
        
        fprintf(ucd_f, "  outer loop\n");
        fprintf(ucd_f, "    vertex %16.10e %16.10e %16.10e\n",finalmesh->node[n0].vert[0],finalmesh->node[n0].vert[1],finalmesh->node[n0].vert[2]);
        fprintf(ucd_f, "    vertex %16.10e %16.10e %16.10e\n",finalmesh->node[n1].vert[0],finalmesh->node[n1].vert[1],finalmesh->node[n1].vert[2]);
        fprintf(ucd_f, "    vertex %16.10e %16.10e %16.10e\n",finalmesh->node[n2].vert[0],finalmesh->node[n2].vert[1],finalmesh->node[n2].vert[2]);
        fprintf(ucd_f, "  endloop\n");
        fprintf(ucd_f, "endfacet\n");
        }
      }
      
    fprintf(ucd_f, "endsolid for_pointwise\n");
    fflush(ucd_f);
    
    fprintf(out_f, "Area summation (x, y, z) = %16.10e %16.10e %16.10e\n",areasum[0],areasum[1],areasum[2]);
    fflush(out_f);
    }
    
  //close ucd file
  fclose(ucd_f);
  
  time(&tm);
  t_char = ctime(&tm);
  fprintf(out_f,"\nFinished with output UCD/STL file at %s.\n",t_char);
  fflush(out_f);

  nonviscbd.destruct();
  
  //***

  fprintf(out_f,"\nWriting restart file.\n");
  fflush(out_f);

  //set up for restart
  //reset file name before restart mesh output
  sname[0] = '\0';
  
  sprintf(sname,"FASTAR_restart.cgns");

  //output vol mesh
  finalmesh->smooth_io(1,geom,sname);
 
  fprintf(out_f,"\nFinished writing restart file.\n");
  fflush(out_f);
  }

  //if (restart == 1 && gentet == 0)
  if (restart == 1)
    {
    fprintf(out_f,"\nReading restart file.\n");
    fflush(out_f);

    //read in file
    //reset file name before volume mesh output
    sname[0] = '\0';
  
    sprintf(sname,"FASTAR_restart.cgns");

    //input vol mesh
    finalmesh->smooth_io(-1,geom,sname);

    fprintf(out_f,"\nFinished reading restart file.\n");
    fflush(out_f); 
    
    printf("\nEnter internal boundary tolerance ->");
    fgets(buff,bdim2,stdin);
    sscanf(buff,"%lg",&tol);
    
    printf("\nEnter geometry boundary tolerance ->");
    fgets(buff,bdim2,stdin);
    sscanf(buff,"%lg",&gtol);
    }

  /*tetgenio *in, *out;
  in = new tetgenio();
  in->initialize();
  out = new tetgenio();
  out->initialize();
  tetgenio::facet *f;
  tetgenio::polygon *p;
  tetgenbehavior flags;

  if (restart == 1 && gentet == 1)
    {
    fprintf(out_f,"\nReading restart file.\n");
    fflush(out_f);

    //read in file
    //reset file name before volume mesh output
    sname[0] = '\0';
  
    sprintf(sname,"FASTAR_restart.cgns");

    //input vol mesh
    finalmesh->smooth_io(-1,geom,sname);

    fprintf(out_f,"\nFinished reading restart file.\n");
    fflush(out_f); 
    
    printf("\nEnter internal boundary tolerance ->");
    fgets(buff,bdim2,stdin);
    sscanf(buff,"%lg",&tol);
    
    printf("\nEnter geometry boundary tolerance ->");
    fgets(buff,bdim2,stdin);
    sscanf(buff,"%lg",&gtol);

    //the main difference here is that I need a map!
    map = new int[finalmesh->nn];

    //init map
    for (i = 0; i < finalmesh->nn; i++)
      map[i] = -1;

    int nni = 0;

    for (i = 0; i < geom->ngb+1; i++)
      {
      for (j = 0; j < finalmesh->boundary[i].elist->max; j++)
        {
        n0 = finalmesh->element[finalmesh->boundary[i].elist->list[j]].f_n[0]->list[0];
        n1 = finalmesh->element[finalmesh->boundary[i].elist->list[j]].f_n[0]->list[1];
        n2 = finalmesh->element[finalmesh->boundary[i].elist->list[j]].f_n[0]->list[2];

        if (map[n0] == -1)
          map[n0] = nni++; 

        if (map[n1] == -1)
          map[n1] = nni++;

        if (map[n2] == -1)
          map[n2] = nni++;
        }
      }

    //set up num facets
    int tempnt = finalmesh->boundary[geom->ngb].elist->max;
    int ngt = 0;

    for (i = 0; i < geom->ngb; i++)
      ngt += finalmesh->boundary[i].elist->max;

    in->firstnumber = 1;
    in->numberofpoints = nni;

    //printf("\npoints=%d\n",in->numberofpoints);

    in->pointlist = new REAL[in->numberofpoints * 3];
    for (i = 0; i < finalmesh->nn; i++)
      if (map[i] >= 0)
        {
        j = map[i]*3;
        in->pointlist[j] = (REAL) finalmesh->node[i].vert[0];
        in->pointlist[j+1] = (REAL) finalmesh->node[i].vert[1];
        in->pointlist[j+2] = (REAL) finalmesh->node[i].vert[2];
        }

    in->numberoffacets = tempnt+ngt;

    //printf("\nfacets=%d\n",in->numberoffacets);    

    in->facetlist = new tetgenio::facet[in->numberoffacets];
    in->facetmarkerlist = new int[in->numberoffacets];
    k = 0; //for numbering
    for (i = 0; i < geom->ngb+1; i++)
      {
      for (j = 0; j < finalmesh->boundary[i].elist->max; j++)
        {
        n0 = finalmesh->element[finalmesh->boundary[i].elist->list[j]].f_n[0]->list[0];
        n1 = finalmesh->element[finalmesh->boundary[i].elist->list[j]].f_n[0]->list[1];
        n2 = finalmesh->element[finalmesh->boundary[i].elist->list[j]].f_n[0]->list[2];
      
        f = &(in->facetlist[k]);
        in->init(f);
        f->numberofpolygons = 1;
        f->polygonlist = new tetgenio::polygon[1];
        f->numberofholes = 0;
        f->holelist = NULL;
        p = &(f->polygonlist[0]);
        in->init(p);
        p->numberofvertices = 3;
        p->vertexlist = new int[p->numberofvertices];
        p->vertexlist[0] = map[n0]+1;
        p->vertexlist[1] = map[n1]+1;
        p->vertexlist[2] = map[n2]+1;
        in->facetmarkerlist[k] = i;
        k++;
        }
      }
    in->save_nodes("FASTAR");
    in->save_poly("FASTAR");

    flags.plc = 1;
    flags.quality = 1;
    flags.quiet = 0;
    flags.epsilon = 1.0e-08;
    flags.nobisect = 2;

    //now, I can delete map
    delete [] map;    

    tetrahedralize(&flags, in, out);

    out->save_nodes("FASTARout");
    out->save_elements("FASTARout");
    out->save_faces("FASTARout");
    
    in->deinitialize(); 
    delete in;
    }

  if (gentet == 1 && restart == 0)
    {
    in->firstnumber = 1;
    in->numberofpoints = nni;

    //printf("\npoints=%d\n",in->numberofpoints);

    in->pointlist = new REAL[in->numberofpoints * 3];
    for (i = 0; i < finalmesh->nn; i++)
      if (map[i] >= 0)
        {
        j = map[i]*3;
        in->pointlist[j] = (REAL) finalmesh->node[i].vert[0];
        in->pointlist[j+1] = (REAL) finalmesh->node[i].vert[1];
        in->pointlist[j+2] = (REAL) finalmesh->node[i].vert[2];
        }

    in->numberoffacets = tempnt+ngt;

    //printf("\nfacets=%d\n",in->numberoffacets);    

    in->facetlist = new tetgenio::facet[in->numberoffacets];
    in->facetmarkerlist = new int[in->numberoffacets];
    k = 0; //for numbering
    for (i = 0; i < geom->ngb+1; i++)
      {
      for (j = 0; j < finalmesh->boundary[i].elist->max; j++)
        {
        n0 = finalmesh->element[finalmesh->boundary[i].elist->list[j]].f_n[0]->list[0];
        n1 = finalmesh->element[finalmesh->boundary[i].elist->list[j]].f_n[0]->list[1];
        n2 = finalmesh->element[finalmesh->boundary[i].elist->list[j]].f_n[0]->list[2];
      
        f = &(in->facetlist[k]);
        in->init(f);
        f->numberofpolygons = 1;
        f->polygonlist = new tetgenio::polygon[1];
        f->numberofholes = 0;
        f->holelist = NULL;
        p = &(f->polygonlist[0]);
        in->init(p);
        p->numberofvertices = 3;
        p->vertexlist = new int[p->numberofvertices];
        p->vertexlist[0] = map[n0]+1;
        p->vertexlist[1] = map[n1]+1;
        p->vertexlist[2] = map[n2]+1;
        in->facetmarkerlist[k] = i;
        k++;
        }
      }
    in->save_nodes("FASTAR");
    in->save_poly("FASTAR");

    flags.plc = 1;
    flags.quality = 1;
    flags.quiet = 0;
    flags.epsilon = 1.0e-08;
    flags.nobisect = 2;

    //now, I can delete map
    delete [] map;    

    tetrahedralize(&flags, in, out);

    out->save_nodes("FASTARout");
    out->save_elements("FASTARout");
    out->save_faces("FASTARout");

    in->deinitialize();
    delete in;
    }*/
    
  //NOTE: Bruce's code will have a full set of bd to be passed in as part of the mesh obj
  //the ngb+1 bd is the internal one
  //he will delete this (as seen later in code) after he has added in tets and nodes, then all we have to do is write out CGNS file from finalmesh
  //we can also dump out UCD creation, but restart seems smart to leave in   

  //NOTICE: smn, smx are smaller or right on as have been upgraded after create_tree
  int num_files = 0, ext = 0, bext = 0, text = 0;
  char **filename;
  
  //vars for tetmesh
  int tetnn = 0, ntet = 0;
  int ntribd = 0;
  int tntribd = 0;
  int **tribd = 0;
  int **tets = 0;
  Point *tetnodes = 0;
  int tempbd = 0, t0, t1, t2, t3, t4;
  double node0, node1, node2;
  
  if (!gentet)
  {
  time(&tm);
  t_char = ctime(&tm);
  fprintf(out_f,"\nBeginning stitching at %s.\n",t_char);
  //now, we need to read in CGNS files from PW
  printf("\nEnter number of tetrahedral mesh files to be stitched to polyhedral mesh ->");
  fgets(buff,bdim2,stdin);
  sscanf(buff,"%d",&num_files);
  
  fprintf(out_f,"\nnum_files to be stitched = %d\n",num_files);
  fflush(out_f);
  
  filename = (char**)malloc(num_files*sizeof(char*));
  for (i = 0; i < num_files; i++)
    {
    filename[i] = (char*)malloc(bdim2*sizeof(char));
    printf("\nEnter filename for file #%d ->",i+1);
    fgets(buff,bdim2,stdin);
    sscanf(buff,"%s",filename[i]);
    fprintf(out_f,"\nStitching Filename #%d = %s\n",i+1,filename[i]);
    fflush(out_f);
    }
  
  //we have to be cognizant that we must append these arrays with info from all files (which allows no duplicate nodes/tris, but may have more nodes)!
  //this is true since will only mesh separate regions when fully distinct
  //not concerned with order, since will all be remapped anyways
  for (i = 0; i < num_files; i++)
    {
    //reset read in
    ext = bext = text = 0;
    
    fprintf(out_f,"\nReading stitching file %d.\n",i+1);
    fflush(out_f);
    
    //open file
    if ((tet_f=fopen(filename[i],"r")) == NULL)
      {
      fprintf(stderr,"\nCouldn't open tet meshfile %s.  Exiting....",filename[i]);
      fflush(stderr);
      exit(0);
      }
      
    //read through all garbage
    for (j = 0; j < 13; j++)
      fgets(buff,bdim2,tet_f);
      
    //now read num nodes
    fgets(buff,bdim2,tet_f);
    sscanf(buff,"%d",&ext);
    
    //now, allocate
    fprintf(out_f,"\nAllocating for %d nodes.\n", tetnn+ext);
    fflush(out_f);
    if (tetnodes == 0)
      tetnodes = (Point*)malloc((tetnn+ext)*sizeof(Point));
    else
      tetnodes = (Point*)realloc((void*)tetnodes, (tetnn+ext)*sizeof(Point));
    
    //read nodes
    for (j = tetnn; j < tetnn+ext; j++)
      {
      fgets(buff,bdim2,tet_f);
      sscanf(buff,"%lg %lg %lg",&(node0),&(node1),&(node2));
      //since not using new, no constructor called, need to read into doubles, init pt
      tetnodes[j] = Point(node0,node1,node2);
      }
      
    //wait to reset tetnn until all remapped
      
    //skip bfaces tag
    fgets(buff,bdim2,tet_f);
    
    //read number, so can loop, but be sure to differentiate before alloc
    fgets(buff,bdim2,tet_f);
    sscanf(buff,"%d",&bext);
    
    tempbd += bext;
    
    fprintf(out_f,"\nFound %d boundary triangles.\n",bext);
    fflush(out_f);
      
    //save temp ntribd
    tntribd = ntribd;
    
    //now, allocate bd arrays
    ntribd += bext;
    
    fprintf(out_f,"\nAllocating for %d boundary triangles.\n",tempbd);
    fflush(out_f);
    
    if (tribd == 0)
      {
      tribd = (int**)calloc(ntribd,sizeof(int*));
      for (k = 0; k < ntribd; k++)
        {
        tribd[k] = (int*)calloc(3,sizeof(int));
        for (l = 0; l < 3; l++)
          tribd[k][l] = -1;
        }
      }
    else
      {
      tribd = (int**)realloc((void*)tribd,ntribd*sizeof(int*));
      for (k = tntribd; k < ntribd; k++)
        {
        tribd[k] = (int*)calloc(3,sizeof(int));
        for (l = 0; l < 3; l++)
          tribd[k][l] = -1;
        }
      }
      
    //now, read in bd
    //set cntr
    k = tntribd; 
    for (j = 0; j < bext; j++)
      {
      fgets(buff,bdim2,tet_f);
      sscanf(buff,"%d %d %d %d %d",&t0,&t1,&t2,&t3,&t4);
      tribd[k][0] = t2-1+tetnn;
      tribd[k][1] = t3-1+tetnn;
      tribd[k][2] = t4-1+tetnn;
      k++;
      }
      
    //finally read in tets
    //skip element line
    fgets(buff,bdim2,tet_f);
    
    //begin counting
    text = 0;
    flag = 0; //hop out when hit vars
    do
      {
      fgets(buff,bdim2,tet_f);
      if (buff[0] == 'V')
        flag = 1;
      else
        text++;
      fgets(buff,bdim2,tet_f); //1's on first line (or Variables), this gets actual tet
      } while (!flag);
      
    fprintf(out_f,"\nAllocating for %d tets.\n",ntet+text);
    fflush(out_f);
    
    if (tets == 0)
      {
      tets = (int**)calloc(text,sizeof(int*));
      for (j = 0; j < text; j++)
        {
        tets[j] = (int*)calloc(4,sizeof(int));
        for (k = 0; k < 4; k++)
          tets[j][k] = -1;
        }
      }
    else
      {
      tets = (int**)realloc((void*)tets,(ntet+text)*sizeof(int*));
      for (j = ntet; j < ntet+text; j++)
        {
        tets[j] = (int*)calloc(4,sizeof(int));
        for (k = 0; k < 4; k++)
          tets[j][k] = -1;
        }
      }
    
    //now, rewind
    rewind(tet_f);
    
    //skip through junk
    for (j = 0; j < 17+ext+bext; j++)
      fgets(buff,bdim2,tet_f);  
    
    //begin reading tets
    for (j = ntet; j < ntet+text; j++)
      {
      fgets(buff,bdim2,tet_f);
      fgets(buff,bdim2,tet_f);
      sscanf(buff,"%d %d %d %d",&(tets[j][0]),&(tets[j][1]),&(tets[j][2]),&(tets[j][3]));
      }
      
    //now, increment tet nodes
    for (j = ntet; j < ntet+text; j++)
      {
      tets[j][0] += tetnn-1;
      tets[j][1] += tetnn-1;
      tets[j][2] += tetnn-1;
      tets[j][3] += tetnn-1;
      }
    
    //now, reset ntet
    ntet += text;
    
    //finally, reset tetnn
    tetnn += ext;
    
    //debug
    //printf("\n%d %d %d %d\n",tets[2][0],tets[2][1],tets[2][2],tets[2][3]);
    //fflush(stdin);
      
    //don't care about vars, so ignore them
    fclose(tet_f);
    
    //now we have all files' info, so we can start stitching!
    }
  
  fprintf(out_f,"\nFinished reading files. Beginning stitching process.\n");
  fflush(out_f); 
  }
  if (gentet)
    {
    time(&tm);
    t_char = ctime(&tm);
    fprintf(out_f,"\nBeginning TETGEN at %s.\n",t_char);
    if (restart == 0)
      delete [] map;
    //need to read in STL file, call tetgen, then read in .nodes and .faces and .elem files, reconstruct mesh
    if (highquality)
      system("./tetgen -pqYYAAT1.0e-12 FASTAR.stl");
    else
      system("./tetgen -pqYYAAT1.0e-12 FASTAR.stl");

    time(&tm);
    t_char = ctime(&tm);
    fprintf(out_f,"\nFinshed TETGEN at %s.\n",t_char);

    //also, create output mesh of only tets
    POLYMESH *tetmesh = new POLYMESH();

    //now, read in files output by tet mesher
    //open node file
    if ((tet_f=fopen("FASTAR.1.node","r")) == NULL)
      {
      fprintf(stderr,"\nCouldn't open tet meshfile FASTAR.1.node.  Exiting....");
      fflush(stderr);
      exit(0);
      }
      
    //read num nodes
    fgets(buff,bdim2,tet_f);
    sscanf(buff,"%d %d %d %d",&tetnn,&i,&i,&i);
    
    //now, allocate
    /*fprintf(out_f,"\nAllocating for %d nodes.\n", out->numberofpoints);
    fflush(out_f);
    tetnodes = (Point*)malloc(out->numberofpoints*sizeof(Point));
    tetmesh->node = new NODE[out->numberofpoints];*/
    fprintf(out_f,"\nAllocating for %d nodes.\n", tetnn);
    fflush(out_f);
    tetnodes = (Point*)malloc(tetnn*sizeof(Point));
    tetmesh->node = new NODE[tetnn];
    tetmesh->nn = tetnn;
    
    //create nodes
    for (i = 0; i < tetnn; i++)
      {
      fgets(buff,bdim2,tet_f);
      sscanf(buff,"%d %lg %lg %lg",&j,&(node0),&(node1),&(node2));
      tetnodes[i] = Point(node0,node1,node2);
      tetmesh->node[i].vert = Point(node0,node1,node2);
      }
    /*for (i = 0; i < out->numberofpoints; i++)
      {
      tetnodes[i] = Point(out->pointlist[3*i],out->pointlist[3*i+1],out->pointlist[3*i+2]);
      tetmesh->node[i].vert = Point(out->pointlist[3*i],out->pointlist[3*i+1],out->pointlist[3*i+2]);
      }

    tetnn = out->numberofpoints;
    tetmesh->nn = out->numberofpoints;*/

    fclose(tet_f);
    tet_f = NULL;

    //open elem file just to get number of tets
    //open face file
    if ((tet_f=fopen("FASTAR.1.ele","r")) == NULL)
      {
      fprintf(stderr,"\nCouldn't open tet meshfile FASTAR.1.ele.  Exiting....");
      fflush(stderr);
      exit(0);
      }
      
    //read num tets
    fgets(buff,bdim2,tet_f);
    sscanf(buff,"%d %d %d",&ntet,&i,&i);
    
    fprintf(out_f,"\nNumber of tets = %d.\n", ntet);
    fflush(out_f);

    fclose(tet_f);
    tet_f = NULL;

    //open face file
    if ((tet_f=fopen("FASTAR.1.face","r")) == NULL)
      {
      fprintf(stderr,"\nCouldn't open tet meshfile FASTAR.1.face.  Exiting....");
      fflush(stderr);
      exit(0);
      }
      
    //read num faces
    fgets(buff,bdim2,tet_f);
    sscanf(buff,"%d %d",&ntribd,&i);
    
    //fprintf(out_f,"\nAllocating for %d facets.\n", out->numberoffacets);
    fprintf(out_f,"\nAllocating for %d facets.\n", ntribd);
    fflush(out_f);

    //tribd = (int**)calloc(out->numberoffacets,sizeof(int*));
    //for (k = 0; k < out->numberoffacets; k++)
    tribd = (int**)calloc(ntribd,sizeof(int*));
    for (k = 0; k < ntribd; k++)
      {
      tribd[k] = (int*)calloc(3,sizeof(int));
      for (l = 0; l < 3; l++)
        tribd[k][l] = -1;
       }

    //ntribd = out->numberoffacets;

    //tetmesh->element = new POLY_ELEMENT[out->numberoffacets + out->numberoftetrahedra];
    tetmesh->element = new POLY_ELEMENT[ntribd + ntet];
    //reset nelems
    tetmesh->nelems = 0;

    //set up bd mesh
    //set number of boundaries...make last one internal bd, go back and delete later
    /*tetmesh->nb = geom->ngb+1; 
  
    //need to add boundaries
    tetmesh->boundary = new BOUNDARY_MESH[geom->ngb + 1];
    
    //first, set all names
    for (i = 0; i < geom->ngb; i++)
      {
      tetmesh->boundary[i].name = new char[400];
      sprintf(finalmesh->boundary[i].name,"%s",geom->g_bname[i+1]);
      }
    
    tetmesh->boundary[geom->ngb].name = new char[400];
    sprintf(finalmesh->boundary[geom->ngb].name,"Stitching Boundary");
     
    //now, set up lists!
    for (i = 0; i < geom->ngb+1; i++)
      {
      tetmesh->boundary[i].elist = new List();
      }*/

    tetmesh->nb = 1; 
  
    //need to add boundaries
    tetmesh->boundary = new BOUNDARY_MESH[1];
    
    tetmesh->boundary[0].name = new char[400];
    sprintf(tetmesh->boundary[0].name,"Generic Boundary");
     
    tetmesh->boundary[0].elist = new List();
 
    //for (j = 0; j < out->numberoffacets; j++)
    for (j = 0; j < ntribd; j++)
      {
      /*f = &(out->facetlist[j]);
      p = &(f->polygonlist[0]);
      tribd[j][0] = p->vertexlist[0] - 1;
      tribd[j][1] = p->vertexlist[1] - 1;
      tribd[j][2] = p->vertexlist[2] - 1;*/

      //read faces
      fgets(buff,bdim2,tet_f);
      sscanf(buff,"%d %d %d %d",&i,&(tribd[j][0]),&(tribd[j][1]),&(tribd[j][2]));
      //decrement
      tribd[j][0]--;
      tribd[j][1]--;
      tribd[j][2]--;

      //set up as elements
      tetmesh->element[tetmesh->nelems].nf = 1;
      tetmesh->element[tetmesh->nelems].f_n = new List*[1];
      tetmesh->element[tetmesh->nelems].f_n[0] = new List();

      tetmesh->element[tetmesh->nelems].f_n[0]->Add_To_List(tribd[j][0]);
      tetmesh->element[tetmesh->nelems].f_n[0]->Add_To_List(tribd[j][1]);
      tetmesh->element[tetmesh->nelems].f_n[0]->Add_To_List(tribd[j][2]);

      //set up as boundarymesh
      //tetmesh->boundary[out->facetmarkerlist[j]].elist->Add_To_List(tetmesh->nelems);
      tetmesh->boundary[0].elist->Add_To_List(tetmesh->nelems);
      
      //inc nelems
      tetmesh->nelems++;
      }

    fclose(tet_f);
    tet_f = NULL;

    //open elem file
    if ((tet_f=fopen("FASTAR.1.ele","r")) == NULL)
      {
      fprintf(stderr,"\nCouldn't open tet meshfile FASTAR.1.ele.  Exiting....");
      fflush(stderr);
      exit(0);
      }
      
    //read num tets
    fgets(buff,bdim2,tet_f);
    sscanf(buff,"%d %d %d",&ntet,&i,&i);
    
    //fprintf(out_f,"\nAllocating for %d tets.\n", out->numberoftetrahedra);
    fprintf(out_f,"\nAllocating for %d tets.\n", ntet);
    fflush(out_f);
  
    //tets = (int**)calloc(out->numberoftetrahedra,sizeof(int*));
    //for (j = 0; j < out->numberoftetrahedra; j++)
    tets = (int**)calloc(ntet,sizeof(int*));
    for (j = 0; j < ntet; j++)
      {
      tets[j] = (int*)calloc(4,sizeof(int));
      for (k = 0; k < 4; k++)
        tets[j][k] = -1;
      } 
    
    //begin reading tets...need only those in region
    int t1, t2, t3, t4, reg, newntet;

    //keeps track of used tets
    newntet = 0;

    //make list of regions
    List reglist;
    reglist.construct();

    for (j = 0; j < numregions; j++)
      reglist.Add_To_List(regions[j]);

    //for (j = 0; j < out->numberoftetrahedra; j++)
    for (j = 0; j < ntet; j++)
      {
      /*tets[j][0] = out->tetrahedronlist[4*j] - 1;
      tets[j][1] = out->tetrahedronlist[4*j+1] - 1;
      tets[j][2] = out->tetrahedronlist[4*j+2] - 1;
      tets[j][3] = out->tetrahedronlist[4*j+3] - 1;*/

      //read faces
      fgets(buff,bdim2,tet_f);
      sscanf(buff,"%d %d %d %d %d %d",&i,&(t1),&(t2),&(t3),&(t4),&(reg));

      if (!reglist.Is_In_List(reg))
        continue;

      tets[newntet][0] = t1-1;
      tets[newntet][1] = t2-1;
      tets[newntet][2] = t3-1;
      tets[newntet][3] = t4-1;

      //set up as elements
      tetmesh->element[tetmesh->nelems].nf = 4;
      tetmesh->element[tetmesh->nelems].f_n = new List*[4];
      for (k = 0; k < 4; k++)
        tetmesh->element[tetmesh->nelems].f_n[k] = new List();

      tetmesh->element[tetmesh->nelems].f_n[0]->Add_To_List(tets[newntet][0]);
      tetmesh->element[tetmesh->nelems].f_n[0]->Add_To_List(tets[newntet][2]);
      tetmesh->element[tetmesh->nelems].f_n[0]->Add_To_List(tets[newntet][1]);

      tetmesh->element[tetmesh->nelems].f_n[1]->Add_To_List(tets[newntet][0]);
      tetmesh->element[tetmesh->nelems].f_n[1]->Add_To_List(tets[newntet][1]);
      tetmesh->element[tetmesh->nelems].f_n[1]->Add_To_List(tets[newntet][3]);

      tetmesh->element[tetmesh->nelems].f_n[2]->Add_To_List(tets[newntet][1]);
      tetmesh->element[tetmesh->nelems].f_n[2]->Add_To_List(tets[newntet][2]);
      tetmesh->element[tetmesh->nelems].f_n[2]->Add_To_List(tets[newntet][3]);

      tetmesh->element[tetmesh->nelems].f_n[3]->Add_To_List(tets[newntet][2]);
      tetmesh->element[tetmesh->nelems].f_n[3]->Add_To_List(tets[newntet][0]);
      tetmesh->element[tetmesh->nelems].f_n[3]->Add_To_List(tets[newntet][3]);
      
      //inc nelems
      newntet++;
      tetmesh->nelems++;
      }

    reglist.destruct();
    //finally, reset ntet
    ntet = newntet;

    fclose(tet_f);
    tet_f = NULL;

    //finally, dump extra nodes and realloc tets
    tets = (int**)realloc((void*)tets,ntet*sizeof(int*));

    //set up map array
    map = new int[tetnn];

    //init
    for (i = 0; i < tetnn; i++)
      map[i] = -1;

    for (i = 0; i < ntet; i++)
     for (j = 0; j < 4; j++)
       map[tets[i][j]] = 0;

    for (i = 0; i < ntribd; i++)
     for (j = 0; j < 3; j++)
       map[tribd[i][j]] = 0;

    j = 0;
    for (i = 0; i < tetnn; i++)
      if (map[i] == 0)
        map[i] = j++;

    //remap tets and bd faces
    for (i = 0; i < ntet; i++)
      for (j = 0; j < 4; j++)
        tets[i][j] = map[tets[i][j]];
 
    for (i = 0; i < ntribd; i++)
      for (j = 0; j < 3; j++)
        tribd[i][j] = map[tribd[i][j]];

    j = 0;
    //finally, reset nodes
    for (i = 0; i < tetnn; i++)
      if (map[i] >= 0)
        {
        tetnodes[map[i]] = Point(tetnodes[i][0],tetnodes[i][1],tetnodes[i][2]);
        j++;
        }

    //reset tetnn
    tetnn = j;

    //realloc
    tetnodes = (Point*)realloc((void*)tetnodes,tetnn*sizeof(Point));

    //now, reset ntet
    //ntet = out->numberoftetrahedra;

    //out->deinitialize();
    //delete out;

    //finally write out and delete tetmesh
    //reset file name before volume mesh output
    sname[0] = '\0';
  
    sprintf(sname,"FASTAR_tet.cgns");

    //output vol mesh
    tetmesh->smooth_io(1,geom,sname);
  
    fprintf(out_f,"\nFinished with output tet CGNS file.\n");
    fflush(out_f);

    delete [] map;

    delete tetmesh;
    } 
  
  //create map of tetnodes...will make elem creation possible, and takes less mem than lists
  int *cmap = new int[tetnn]; 
  //init to -1
  for (i = 0; i < tetnn; i++)
    cmap[i] = -1;
   
  fprintf(out_f,"\nSetting up octree sorted bd facet list.\n");
  fflush(out_f);  
      
  //will need to set up all in same octree
  //create Octree sorted facet list
  Octree_Storage *Octree_root;
  POLY_ELEMENT *pntr = 0;
  PList *elist = 0;
  double stol[3]; //will hold slop factor always
  double pt[3]; //will hold pt for octree retrieve
  int nd0, nd1, nd2, flag1 = 0;
  Point p2;

  List nonviscbd;
  nonviscbd.construct();

  for (i = 0; i < finalmesh->nb && viscous; i++)
     {
     flag = 0; 
     for (j = 0; j < nvbd && !flag; j++)
       {
       if (i+1 == vbd[j])
         flag++;
       }
     if (!flag)
       {
       nonviscbd.Add_To_List(i);
       }
     }

  nonviscbd.print(out_f);
  
  //reset points
  p0 = Point(1.0e20, 1.0e20, 1.0e20);
  p1 = Point(-1.0e20, -1.0e20, -1.0e20);
  
  int n_facets = 0; //to store total number
  
  //find extents for geom
  for (i=0; i < finalmesh->nb; i++)
    {
    if (!nonviscbd.Is_In_List(i) && viscous)
      continue;
    for (j=0; j < finalmesh->boundary[i].elist->max; j++)
      {
      n0 = finalmesh->element[finalmesh->boundary[i].elist->list[j]].f_n[0]->list[0];
      n1 = finalmesh->element[finalmesh->boundary[i].elist->list[j]].f_n[0]->list[1];
      n2 = finalmesh->element[finalmesh->boundary[i].elist->list[j]].f_n[0]->list[2];
      
      p0[0] = MIN(finalmesh->node[n2].vert[0],MIN(finalmesh->node[n1].vert[0],MIN(p0[0],finalmesh->node[n0].vert[0])));
      p0[1] = MIN(finalmesh->node[n2].vert[1],MIN(finalmesh->node[n1].vert[1],MIN(p0[1],finalmesh->node[n0].vert[1])));
      p0[2] = MIN(finalmesh->node[n2].vert[2],MIN(finalmesh->node[n1].vert[2],MIN(p0[2],finalmesh->node[n0].vert[2])));
      p1[0] = MAX(finalmesh->node[n2].vert[0],MAX(finalmesh->node[n1].vert[0],MAX(p1[0],finalmesh->node[n0].vert[0])));
      p1[1] = MAX(finalmesh->node[n2].vert[1],MAX(finalmesh->node[n1].vert[1],MAX(p1[1],finalmesh->node[n0].vert[1])));
      p1[2] = MAX(finalmesh->node[n2].vert[2],MAX(finalmesh->node[n1].vert[2],MAX(p1[2],finalmesh->node[n0].vert[2])));
      n_facets++;
      }
    }
    
  for (j = 0; j < 3; j++)
    {
    olo[j] = p0[j];
    ohi[j] = p1[j];
    }
    
  fprintf(out_f,"\nFound %d bd facet extents.\n",n_facets);
  fflush(out_f);

  //set for mom
  Octree_Storage *dummy=0;
  //set up tlo, thi
  double (*tlo)[3], (*thi)[3];
  //set to reuse and pass in
  void* *ptr;
  //create storage
  Octree_root = new Octree_Storage((Octree_Storage*)dummy,olo,ohi);
  
  //set up ptrs
  ptr = new void*[n_facets];
  tlo = new double[n_facets][3];
  thi = new double[n_facets][3];
  
  k = 0; //set for facet numbering
  //store in octree
  for (i=0; i < finalmesh->nb; i++)
    {
    if (!nonviscbd.Is_In_List(i) && viscous)
      continue;
    for (j=0; j < finalmesh->boundary[i].elist->max; j++)
      {
      n0 = finalmesh->element[finalmesh->boundary[i].elist->list[j]].f_n[0]->list[0];
      n1 = finalmesh->element[finalmesh->boundary[i].elist->list[j]].f_n[0]->list[1];
      n2 = finalmesh->element[finalmesh->boundary[i].elist->list[j]].f_n[0]->list[2];
      
      tlo[k][0] = MIN(finalmesh->node[n2].vert[0],MIN(finalmesh->node[n1].vert[0],finalmesh->node[n0].vert[0]));
      tlo[k][1] = MIN(finalmesh->node[n2].vert[1],MIN(finalmesh->node[n1].vert[1],finalmesh->node[n0].vert[1]));
      tlo[k][2] = MIN(finalmesh->node[n2].vert[2],MIN(finalmesh->node[n1].vert[2],finalmesh->node[n0].vert[2]));
      thi[k][0] = MAX(finalmesh->node[n2].vert[0],MAX(finalmesh->node[n1].vert[0],finalmesh->node[n0].vert[0]));
      thi[k][1] = MAX(finalmesh->node[n2].vert[1],MAX(finalmesh->node[n1].vert[1],finalmesh->node[n0].vert[1]));
      thi[k][2] = MAX(finalmesh->node[n2].vert[2],MAX(finalmesh->node[n1].vert[2],finalmesh->node[n0].vert[2]));
      ptr[k] = (void*)&(finalmesh->element[finalmesh->boundary[i].elist->list[j]]);
      k++;
      }
    }
    
  if (k != n_facets)
    {
    fprintf(out_f, "\nYou have %d facets, but you only processed %d into Octree.  Exiting....\n",n_facets,k);
    fflush(out_f);
    exit(0);
    }
  
  //actually store facets
  Octree_root->Store_In_Octree(n_facets,ptr,tlo,thi);
  
  fprintf(out_f, "\nDone storing facets in Octree.\n");
  fflush(out_f);
     
  //free mem
  delete[] ptr;
  delete[] tlo;
  delete[] thi;
  nonviscbd.destruct();
  
  //be sure pntr reset
  pntr = 0;
  
  //create new Plist
  elist = new PList();
  //dimension a priori, reset max within
  elist->Redimension(n_facets);
    
  fprintf(out_f,"\nRemapping boundary nodes using Octree, so can remap tets and add nodes.\n");
  fflush(out_f);
  
  //set tolarance as lower of geom and inner tol
  tol = MIN(tol,gtol)*1.0e-02;
 
  //find using octree for remapping nodes for tet creation ONLY..all geom in mesh obj and will dump internal
  for (i = 0; i < ntribd; i++)
    {
    n0 = tribd[i][0];
    n1 = tribd[i][1];
    n2 = tribd[i][2];
    
    //find centroid
    cent = (tetnodes[n0] + tetnodes[n1] + tetnodes[n2])/3.0;
    
    //set level to -1 so returns all (only reason to go thru levels is to change slop factor) 
    l = -1;
    
    //now, redim elist
    elist->max = 0;
    
    //set up pt, reset slop
    for (k = 0; k < 3; k++)
      {
      pt[k] = cent[k]; 
      stol[k] = -1.0e+20;
      }
    
    //set up slop factor
    for (k = 0; k < 3; k++)
      {
      stol[k] = MIN(MAX(stol[k],MAX(0.5*fabs(tetnodes[n1][k] - tetnodes[n0][k]),ctol)),1.0e-02);
      stol[k] = MIN(MAX(stol[k],MAX(0.5*fabs(tetnodes[n2][k] - tetnodes[n1][k]),ctol)),1.0e-02);
      stol[k] = MIN(MAX(stol[k],MAX(0.5*fabs(tetnodes[n0][k] - tetnodes[n2][k]),ctol)),1.0e-02);
      }
    
    //now, receive list of facets in the neighborhood to cut
    Octree_root->retrieve_list(pt,stol,l,elist);
    
    if (elist->max == 0)
      {
      fprintf(stderr, "\nYou failed to find candidate boundary triangles for tri %d from the original list.  Exiting....\n",i);
      fflush(stderr);
      exit(0);
      }
    
    flag = 0; //jump out once found, check that is found
    for (k = 0; k < elist->max && !flag; k++)
      {
      flag1 = 0;
      //set pntr pulled out of elist to what it actually is...a facet
      pntr = (POLY_ELEMENT*)elist->list[k];
      
      //set up nodes from facet
      nd0 = pntr->f_n[0]->list[0];
      nd1 = pntr->f_n[0]->list[1];
      nd2 = pntr->f_n[0]->list[2];
      p0 = finalmesh->node[nd0].vert;
      p1 = finalmesh->node[nd1].vert;
      p2 = finalmesh->node[nd2].vert;
      
      //check for correct nodes, set cmap...if not all found, will kick out
      if (distance(p0, tetnodes[n0]) < tol)
        {
        cmap[n0] = nd0;
        flag1++;
        }
      else if (distance(p0, tetnodes[n1]) < tol)
        {
        cmap[n1] = nd0;
        flag1++;
        }
      else if (distance(p0, tetnodes[n2]) < tol)
        {
        cmap[n2] = nd0;
        flag1++;
        }
        
      if (distance(p1, tetnodes[n0]) < tol)
        {
        cmap[n0] = nd1;
        flag1++;
        }
      else if (distance(p1, tetnodes[n1]) < tol)
        {
        cmap[n1] = nd1;
        flag1++;
        }
      else if (distance(p1, tetnodes[n2]) < tol)
        {
        cmap[n2] = nd1;
        flag1++;
        }
        
      if (distance(p2, tetnodes[n0]) < tol)
        {
        cmap[n0] = nd2;
        flag1++;
        }
      else if (distance(p2, tetnodes[n1]) < tol)
        {
        cmap[n1] = nd2;
        flag1++;
        }
      else if (distance(p2, tetnodes[n2]) < tol)
        {
        cmap[n2] = nd2;
        flag1++;
        }
        
      if (flag1 == 3)
        flag = 1;
      }
      
    //check extant
    if (!flag)
      {
      fprintf(stderr, "\nYou failed to find a boundary triangle for tri %d from the original list. Flag = %d.  Exiting....\n",i,flag1);
      fflush(stderr);
      exit(0);
      }
    }
    
  //finally, free octree mem 
  dummy = 0; 
  if (elist != 0)
    delete elist;
  if (Octree_root != 0)
    delete Octree_root;
  
  fprintf(out_f,"\nDone remapping facets.  Adding in tetmesher generated nodes.\n");
  fflush(out_f);
  
  //finally add in pw generated nodes
  for (i = 0; i < tetnn; i++)
    {
    if (cmap[i] >= 0)
      continue; //not pw node
    //check for alloc
    if ((finalmesh->nn+1) > finalmesh->node_dim)
      my_mem(finalmesh->node_dim, finalmesh->nn, &(finalmesh->node), (finalmesh->nn+1), NODE_CHUNK);
    finalmesh->node[finalmesh->nn].vert = tetnodes[i];
    //set map
    cmap[i] = finalmesh->nn;
    //inc finalmesh->nn
    finalmesh->nn++;
    }
    
  fprintf(out_f,"\nDeleting internal facets from final polymesh.\n");
  fflush(out_f);
    
  //now, we need to dump bd elements from internal bd
  //no need to realloc since there will be at least this many tets getting added in, and these were last elem added in...just change nelem
  fprintf(out_f,"\nElements = %d.\n", finalmesh->nelems);
  fflush(out_f); 

  j = 0;
  for (i = 0; i < finalmesh->nelems; i++)
    if (finalmesh->element[i].nf == 0)
      j++;

  fprintf(out_f,"\nnumber of null elem = %d\n",j);
  fflush(out_f);

   j = 0;
  for (i = 0; i < finalmesh->nelems; i++)
    if (finalmesh->element[i].nf == 1)
      j++;

  fprintf(out_f,"\nnumber of bd elem = %d\n",j);
  fflush(out_f);

  if (viscous)
    {
    for (i = 0; i < finalmesh->boundary[geom->ngb+1].elist->max; i++)
      {
      finalmesh->element[finalmesh->boundary[geom->ngb+1].elist->list[i]].f_n[0]->Redimension(0);
      delete finalmesh->element[finalmesh->boundary[geom->ngb+1].elist->list[i]].f_n[0];
      finalmesh->element[finalmesh->boundary[geom->ngb+1].elist->list[i]].nf = 0;
      }
    }
  
  for (i = 0; i < finalmesh->boundary[geom->ngb].elist->max; i++)
    {
    finalmesh->element[finalmesh->boundary[geom->ngb].elist->list[i]].f_n[0]->Redimension(0);
    delete finalmesh->element[finalmesh->boundary[geom->ngb].elist->list[i]].f_n[0];
    finalmesh->element[finalmesh->boundary[geom->ngb].elist->list[i]].nf = 0;
    }

  j = 0;
  for (i = 0; i < finalmesh->nelems; i++)
    if (finalmesh->element[i].nf == 0)
      j++;

  fprintf(out_f,"\nnumber of null elem = %d\n",j);
  fflush(out_f);

   j = 0;
  for (i = 0; i < finalmesh->nelems; i++)
    if (finalmesh->element[i].nf == 1)
      j++;

  fprintf(out_f,"\nnumber of bd elem = %d\n",j);
  fflush(out_f);
    
  //now, decrement nelems
  finalmesh->nelems -= finalmesh->boundary[geom->ngb].elist->max;
  if (viscous)
    finalmesh->nelems -= finalmesh->boundary[geom->ngb+1].elist->max;
  
  //now, redim that
  finalmesh->boundary[geom->ngb].elist->Redimension(0);
  if (viscous)
    finalmesh->boundary[geom->ngb+1].elist->Redimension(0);
  
  //delete that
  delete finalmesh->boundary[geom->ngb].elist;
  if (viscous)
    delete finalmesh->boundary[geom->ngb+1].elist;
  
  //reset nb
  finalmesh->nb--;
  if (viscous)
    finalmesh->nb--;

  fprintf(out_f,"\nAfter deleting bd elem, elements = %d.\n", finalmesh->nelems);
  fflush(out_f); 
    
  fprintf(out_f,"\nAllocating for new tet elements.\n");
  fflush(out_f);  
    
  //finally, add in tet elements
  //alloc
  if ((finalmesh->nelems+ntet) > finalmesh->element_dim)
    {
    my_mem(finalmesh->element_dim, finalmesh->nelems, &(finalmesh->element), (finalmesh->nelems+ntet), ELEMENT_CHUNK);
    }
  for (k = finalmesh->nelems; k < (finalmesh->nelems+ntet); k++)
    finalmesh->element[k].Initialize();
   
  fprintf(out_f,"\nCreating/mapping new tet elements.\n");
  fflush(out_f);        
                        
  //now, add in each tet, which will always have four faces
  //THIS ASSUMES TETS HAVE BEEN REWOUND BY CONVERT TO PROPER ORIENTATION
  for (i = 0; i < ntet; i++)
    {
    //add in face lists
    finalmesh->element[finalmesh->nelems].nf = 4;
    finalmesh->element[finalmesh->nelems].f_n = new List*[4];
    for (j = 0; j < 4; j++)
      finalmesh->element[finalmesh->nelems].f_n[j] = new List();
      
    //create each face in CGNS order
    for (j = 0; j < 4; j++)
      {
      switch (j) 
        {
        case 0:
          n0 = tets[i][0];
          n1 = tets[i][2];
          n2 = tets[i][1];
        break;
        case 1:
          n0 = tets[i][0];
          n1 = tets[i][1];
          n2 = tets[i][3];
        break;
        case 2:
          n0 = tets[i][1];
          n1 = tets[i][2];
          n2 = tets[i][3];
        break;
        case 3:
          n0 = tets[i][2];
          n1 = tets[i][0];
          n2 = tets[i][3];
        break;
        default:
        break;
        }
    
      finalmesh->element[finalmesh->nelems].f_n[j]->Add_To_List(cmap[n0]);
      finalmesh->element[finalmesh->nelems].f_n[j]->Add_To_List(cmap[n1]);
      finalmesh->element[finalmesh->nelems].f_n[j]->Add_To_List(cmap[n2]);
      }
    
    finalmesh->nelems++;
    }

  fprintf(out_f,"\n#Elem = %d\n",finalmesh->nelems);
  fflush(out_f);

  j = 0;
  for (i = 0; i < finalmesh->nelems; i++)
    if (finalmesh->element[i].nf == 4)
      j++;

  fprintf(out_f,"\nnumber of tet elem = %d\n",j);
  fflush(out_f);

  j = 0;
  for (i = 0; i < finalmesh->nelems; i++)
    if (finalmesh->element[i].nf == 1)
      j++;

  fprintf(out_f,"\nnumber of bd elem = %d\n",j);
  fflush(out_f);

  j = 0;
  for (i = 0; i < finalmesh->nelems; i++)
    if (finalmesh->element[i].nf == 5)
      j++;

  fprintf(out_f,"\nnumber of prism elem = %d\n",j);
  fflush(out_f);

  j = 0;
  for (i = 0; i < finalmesh->nelems; i++)
    if (finalmesh->element[i].nf == 0)
      j++;

  fprintf(out_f,"\nnumber of null elem = %d\n",j);
  fflush(out_f);

  j = 0;
  for (i = 0; i < finalmesh->nelems; i++)
    if (finalmesh->element[i].nf > 6)
      j++;

  fprintf(out_f,"\nnumber of poly elem = %d\n",j);
  fflush(out_f);
    
  //no longer need cmap
  delete [] cmap;
  
  time(&tm);
  t_char = ctime(&tm);
  fprintf(out_f,"\nFinished stitching process.  About to write out final CGNS file at %s.\n",t_char);
  fflush(out_f);
  
  //finally, delete stitching mem
  free(tetnodes);
  for (i = 0; i < ntet; i++)
    free(tets[i]);
  free(tets);
  for (j = 0; j < ntribd; j++)
    {
    free(tribd[j]);
    }
  free(tribd);
  
#endif
  
  //reset file name before volume mesh output
  sname[0] = '\0';
  
  sprintf(sname,"FASTAR.cgns");

  //output vol mesh
  finalmesh->smooth_io(1,geom,sname);
  
  fprintf(out_f,"\nFinished with output CGNS file.\n");
  fflush(out_f);

  //eventually will add in dynamic here.  Then, fix early otmesh deletion, delete finalmesh, start movement process, keeping delete otmesh at the very end of the process
  //would need calls to P_VLI, ux, spacing library first, then, delete finalmesh, then move mesh with trig, find needed spots, regenerate needed spots, begin again...or perhaps use P_SMOOTH to move instead!
  //NOTICE: smn, smx are smaller or right on as have been upgraded after create_tree
  //delete otmesh;
  //NOTE: due to my_mem calls, destructor has been deprecated and should not be trusted!
  delete finalmesh;
  delete ovst;
      
  return;
}
  
void POLYMESH::set_up_tri(int &tnt, List *tnode, int *map)
  {
  //set up vars
  Vector dr, norm, v1, v2;
  norm = Vector(0.0,0.0,0.0);
  int n;
  
  //debug
  //for (n = 0; n < tnode->max; n++)
    //if (tnode->list[n] == -1)
      //fprintf(out_f, "\n-1 in tnode tnt = %d\n",tnt);
  //fflush(out_f);
            
  //find cg of face
  Point cg = Point(0.0,0.0,0.0);
  for (n = 0; n < tnode->max; n++)
    cg += node[tnode->list[n]].vert;
  cg /= tnode->max;
  //set up trimesh arrays
  int *nbs, ***bs;
  nbs = new int[1];
  bs = new int**[1];
  bs[0] = new int*[tnode->max];
  nbs[0] = 0;
  for (n = 0; n < tnode->max-1; n++)
    {
    bs[0][nbs[0]] = new int[2];
    bs[0][nbs[0]][0] = tnode->Index(tnode->list[n]);
    bs[0][nbs[0]][1] = tnode->Index(tnode->list[n+1]);
    nbs[0]++;
    v1 = Vector(cg,node[tnode->list[n]].vert);
    v2 = Vector(cg,node[tnode->list[n+1]].vert);
    norm += v1 % v2;
    }
  //catch connector edge
  bs[0][nbs[0]] = new int[2];
  bs[0][nbs[0]][0] = tnode->Index(tnode->list[tnode->max-1]);
  bs[0][nbs[0]][1] = tnode->Index(tnode->list[0]);
  nbs[0]++;
  v1 = Vector(cg,node[tnode->list[tnode->max-1]].vert);
  v2 = Vector(cg,node[tnode->list[0]].vert);
  norm += v1 % v2;
  //finally, normalize
  norm.normalize();
  
  double *x, *y;
  x = new double[tnode->max];
  y = new double[tnode->max];
  
  //get list0 from above..this allows us to set up x and y w/out respect to which plane is constant 
  v1 = Vector(cg,node[tnode->list[0]].vert);
  double dot = v1 * norm;
  v1 -= norm*dot;
  v1.normalize();
  v2 = norm % v1;
  v2.normalize();
  double dsmin = 1.0e20;
  for (n = 0; n < tnode->max-1; n++)
    {
    dr = Vector(cg,node[tnode->list[n]].vert);
    int i0 = tnode->Index(tnode->list[n]);
    x[i0] = dr*v1;
    y[i0] = dr*v2;
    dr = Vector(cg,node[tnode->list[n+1]].vert);
    int i1 = tnode->Index(tnode->list[n+1]);
    x[i1] = dr*v1;
    y[i1] = dr*v2;
    double dx = x[i1]-x[i0];
    double dy = y[i1]-y[i0];
    double mag = sqrt(dx*dx+dy*dy);
    dsmin = MIN(dsmin,mag);
    }
  //do end around
  dr = Vector(cg,node[tnode->list[tnode->max-1]].vert);
  int i0 = tnode->max-1;
  x[i0] = dr*v1;
  y[i0] = dr*v2;
  dr = Vector(cg,node[tnode->list[0]].vert);
  int i1 = 0;
  x[i1] = dr*v1;
  y[i1] = dr*v2;
  double dx = x[i1]-x[i0];
  double dy = y[i1]-y[i0];
  double mag = sqrt(dx*dx+dy*dy);
  dsmin = MIN(dsmin,mag);
  
  /*fprintf(out_f, "\n\ntnt = %d\n",tnt);
  for (n = 0; n < tnode->max; n++)
    node[tnode->list[n]].vert.print(out_f);
  fflush(out_f);*/
  
  if (dsmin < 1.0e-15)  // coincident points created in projection, use distribution around unit circle
    {
    fprintf(out_f,"\nset_up_tri: coincident points detected in triangulation process!");
    fflush(out_f);
    //MPI_Abort(MPI_COMM_WORLD,0);
    exit(0);
    }
  
  int tdim = tnode->max*5;
  int (*trit)[3] = new int[tdim][3];
  int lnt, t;
  //call trimesh
  
  if ((lnt = trimesh(tnode->max,tdim,1,0,nbs,bs,x,y,trit)) > 0)
    {
    tnt+= lnt;
    //now, need to create tris
    if ((tnt+nelems) > element_dim)
      {
      //allocate mem
      my_mem(element_dim, nelems, &(element), tnt+nelems, ELEMENT_CHUNK);
      }
    for (n = nelems; n < (nelems+tnt); n++)
      {
      element[n].Initialize();
      element[n].nf = 1;
      element[n].f_n = new List*[1];
      element[n].f_n[0] = new List();
      }
    //take the resulting triangles and set up tri
    for (t=0; t < lnt; t++)
      {
      for (n = 0; n < 3; n++)
        {
        //set up tri
        element[nelems+t].f_n[0]->Add_To_List(tnode->list[trit[t][n]]);
        if (!overset)
          map[tnode->list[trit[t][n]]] = -2;
        }
      //set bd type
      if (viscous && overset)
        boundary[0].elist->Add_To_List(nelems+t); 
      else if (viscous && !overset)
        boundary[nb - 2].elist->Add_To_List(nelems+t); 
      else
        boundary[nb - 1].elist->Add_To_List(nelems+t);
      }
    } 
  else
    {
    fprintf(stderr,"\nTrimesh failed.\n");
    fflush(stderr);
    //MPI_Abort(MPI_COMM_WORLD,0);
    exit(0);
    }
  
  //free mem
  delete[] trit;
  delete[] x;
  delete[] y;
  for (t=0; t < nbs[0]; t++)
    delete[] bs[0][t];
  delete[] nbs;
  delete[] bs[0];
  delete[] bs;
  return;
  }

void POLYMESH::set_up_tri_elem(int &tnt, List *tnode, int *map)
  {
  //set up vars
  Vector dr, norm, v1, v2;
  norm = Vector(0.0,0.0,0.0);
  int n;
  
  //debug
  //for (n = 0; n < tnode->max; n++)
    //if (tnode->list[n] == -1)
      //fprintf(out_f, "\n-1 in tnode tnt = %d\n",tnt);
  //fflush(out_f);
            
  //find cg of face
  Point cg = Point(0.0,0.0,0.0);
  for (n = 0; n < tnode->max; n++)
    cg += node[tnode->list[n]].vert;
  cg /= tnode->max;
  //set up trimesh arrays
  int *nbs, ***bs;
  nbs = new int[1];
  bs = new int**[1];
  bs[0] = new int*[tnode->max];
  nbs[0] = 0;
  for (n = 0; n < tnode->max-1; n++)
    {
    bs[0][nbs[0]] = new int[2];
    bs[0][nbs[0]][0] = tnode->Index(tnode->list[n]);
    bs[0][nbs[0]][1] = tnode->Index(tnode->list[n+1]);
    nbs[0]++;
    v1 = Vector(cg,node[tnode->list[n]].vert);
    v2 = Vector(cg,node[tnode->list[n+1]].vert);
    norm += v1 % v2;
    }
  //catch connector edge
  bs[0][nbs[0]] = new int[2];
  bs[0][nbs[0]][0] = tnode->Index(tnode->list[tnode->max-1]);
  bs[0][nbs[0]][1] = tnode->Index(tnode->list[0]);
  nbs[0]++;
  v1 = Vector(cg,node[tnode->list[tnode->max-1]].vert);
  v2 = Vector(cg,node[tnode->list[0]].vert);
  norm += v1 % v2;
  //finally, normalize
  norm.normalize();
  
  double *x, *y;
  x = new double[tnode->max];
  y = new double[tnode->max];
  
  //get list0 from above..this allows us to set up x and y w/out respect to which plane is constant 
  v1 = Vector(cg,node[tnode->list[0]].vert);
  double dot = v1 * norm;
  v1 -= norm*dot;
  v1.normalize();
  v2 = norm % v1;
  v2.normalize();
  double dsmin = 1.0e20;
  for (n = 0; n < tnode->max-1; n++)
    {
    dr = Vector(cg,node[tnode->list[n]].vert);
    int i0 = tnode->Index(tnode->list[n]);
    x[i0] = dr*v1;
    y[i0] = dr*v2;
    dr = Vector(cg,node[tnode->list[n+1]].vert);
    int i1 = tnode->Index(tnode->list[n+1]);
    x[i1] = dr*v1;
    y[i1] = dr*v2;
    double dx = x[i1]-x[i0];
    double dy = y[i1]-y[i0];
    double mag = sqrt(dx*dx+dy*dy);
    dsmin = MIN(dsmin,mag);
    }
  //do end around
  dr = Vector(cg,node[tnode->list[tnode->max-1]].vert);
  int i0 = tnode->max-1;
  x[i0] = dr*v1;
  y[i0] = dr*v2;
  dr = Vector(cg,node[tnode->list[0]].vert);
  int i1 = 0;
  x[i1] = dr*v1;
  y[i1] = dr*v2;
  double dx = x[i1]-x[i0];
  double dy = y[i1]-y[i0];
  double mag = sqrt(dx*dx+dy*dy);
  dsmin = MIN(dsmin,mag);
  
  /*fprintf(out_f, "\n\ntnt = %d\n",tnt);
  for (n = 0; n < tnode->max; n++)
    node[tnode->list[n]].vert.print(out_f);
  fflush(out_f);*/
  
  if (dsmin < 1.0e-15)  // coincident points created in projection, use distribution around unit circle
    {
    fprintf(out_f,"\nset_up_tri: coincident points detected in triangulation process!");
    fflush(out_f);
    //MPI_Abort(MPI_COMM_WORLD,0);
    exit(0);
    }
  
  int tdim = tnode->max*5;
  int (*trit)[3] = new int[tdim][3];
  int lnt, t;
  //call trimesh
  
  if ((lnt = trimesh(tnode->max,tdim,1,0,nbs,bs,x,y,trit)) > 0)
    {
    tnt+= lnt;
    //now, need to create tris
    if ((tnt+nelems) > element_dim)
      {
      //allocate mem
      my_mem(element_dim, nelems, &(element), tnt+nelems, ELEMENT_CHUNK);
      }
    for (n = nelems; n < (nelems+tnt); n++)
      {
      element[n].Initialize();
      element[n].nf = 1;
      element[n].f_n = new List*[1];
      element[n].f_n[0] = new List();
      }
    //take the resulting triangles and set up tri
    for (t=0; t < lnt; t++)
      {
      for (n = 0; n < 3; n++)
        {
        //set up tri
        element[nelems].f_n[0]->Add_To_List(tnode->list[trit[t][n]]);
        map[tnode->list[trit[t][n]]] = -2;
        }
      //set bd type
      boundary[nb - 1].elist->Add_To_List(nelems); 
      nelems++;
      }
    } 
  else
    {
    fprintf(stderr,"\nTrimesh failed.\n");
    fflush(stderr);
    //MPI_Abort(MPI_COMM_WORLD,0);
    exit(0);
    }
  
  //free mem
  delete[] trit;
  delete[] x;
  delete[] y;
  for (t=0; t < nbs[0]; t++)
    delete[] bs[0][t];
  delete[] nbs;
  delete[] bs[0];
  delete[] bs;
  return;
  }
