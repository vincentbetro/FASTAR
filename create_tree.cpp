#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "Point.h"
#include "Vector.h"
#include "Util.h"
#include "split_tree.h"
#include "Spacing_Field.h"

#define ABS(x) ((x) >= 0 ? (x) : -(x))
#define SIGN(a,b) ((b) < 0.0 ? -ABS(a) : ABS(a))

//global variables
extern int my_rank; /*rank of process*/
extern int num_procs; /*number of processes*/
extern int adapt; //using adaptation or not
extern int ugeom; //spacing file only (1) or both spacing and geom used
extern int evencut; //if we want nice box around geometry
extern int bb; //num bds to be cut
extern int viscous; //insert layersextern
extern double vspace; //viscous spacing
extern int nvbd; //number viscous bd
extern int highquality; //allows desirable constraints to be enforced
extern int altflood; //type of flood fill
extern int floodnum; //cube root of num voxels to ray cast
extern double factor; //geom prog factor
extern FILE *in_f, *jou_f, *out_f, *debug_f; /*input output journal files global*/
double sum_layer(int nl, double geomprog);
//split key
//0 none, 1 x, 2 y, 3 z
//returns nvoxels
void split_tree::refine_voxel(int &new_voxel, int k, int code)
  { 
  int i, j, tdim, tnv;
  
  //should never get here
  if (code == 0)
    return;
  
  //prep for new voxels
  tnv = new_voxel;
  tnv+= 2;

  if (((tnv-1)%100000) == 0)
    {
    fprintf(out_f,"\nCurrent nvoxels = %d\n",tnv);
    fflush(out_f);
    }
  
  //now, need to create next level of voxels if not set up yet
  if (tnv > voxel_dim)
    {
    //set temp voxel_dim
    tdim = voxel_dim;
    //allocate mem
    //fprintf(out_f,"\nvoxel_dim = %d, nvoxels = %d, new_voxel = %d, tnv = %d\n",voxel_dim,nvoxels,new_voxel,tnv);
    my_mem(voxel_dim, nvoxels, &voxels, tnv, STCHUNK);
    //fprintf(out_f,"\nvoxel_dim = %d, nvoxels = %d, new_voxel = %d, tnv = %d\n",voxel_dim,nvoxels,new_voxel,tnv);
    //fflush(out_f);
    //set mom to -1, split/cut to 0
    for (i = tdim; i < voxel_dim; i++)
      {
      voxels[i].split = 0;
      voxels[i].mom = -1;
      voxels[i].cut = 0;
      voxels[i].children.max = 0;
      for (j = 0; j < 27; j++)
        voxels[i].nodes[j] = -1;
      }
    }
    
  //set new_voxel
  new_voxel = tnv;
    
  //set construct since from a class and dim init to random number, causing free/realloc and seg fault
  voxels[k].children.construct();      
  voxels[k].children.Add_To_List(new_voxel-2);
  voxels[k].children.Add_To_List(new_voxel-1);
  
  //set up points for new cornerpts
  Point *pt = new Point[6];
  Point p0, p1;
  double begx, begy, begz, endx, endy, endz, halfx, halfy, halfz;
  
  //calculate physical nodes and place them in array
  p0 = voxels[k].cornerpts[0];
  p1 = voxels[k].cornerpts[1];
    
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
          
  //set up pts of orig voxel
  //8,11,12,14,17,18
  pt[0] = Point(halfx,begy,begz);
  pt[1] = Point(begx,halfy,begz);
  pt[2] = Point(begx,begy,halfz);
  pt[3] = Point(endx,endy,halfz);
  pt[4] = Point(endx,halfy,endz);
  pt[5] = Point(halfx,endy,endz);
  
  //set split
  voxels[k].split = code;
  
  //based on code, decide how to refine
  switch (code) 
    {
      case 1:
          //(cornerpts[0],18,8,cornerpts[1])
          //set up corner points for new voxels, with "left" as 0 , also set mom
          voxels[voxels[k].children.list[0]].cornerpts[0] = p0;
          voxels[voxels[k].children.list[0]].cornerpts[1] = pt[5];
          voxels[voxels[k].children.list[1]].cornerpts[0] = pt[0];
          voxels[voxels[k].children.list[1]].cornerpts[1] = p1;
          break;
      case 2:
          //(cornerpts[0],17,11,cornerpts[1])
          //set up corner points for new voxels, with "front" as 0 , also set mom
          voxels[voxels[k].children.list[0]].cornerpts[0] = p0;
          voxels[voxels[k].children.list[0]].cornerpts[1] = pt[4];
          voxels[voxels[k].children.list[1]].cornerpts[0] = pt[1];
          voxels[voxels[k].children.list[1]].cornerpts[1] = p1;
          break;
      case 3:
          //(cornerpts[0],14,12,cornerpts[1])
          //set up corner points for new voxels, with "bottom" as 0 , also set mom
          voxels[voxels[k].children.list[0]].cornerpts[0] = p0;
          voxels[voxels[k].children.list[0]].cornerpts[1] = pt[3];
          voxels[voxels[k].children.list[1]].cornerpts[0] = pt[2];
          voxels[voxels[k].children.list[1]].cornerpts[1] = p1;
          break;
      default:
          fprintf(stderr,"\nYou have an invalid code.  Exiting....\n");
          fflush(stderr);
          exit(0);
          break;
    }
    
  //set mom
  voxels[voxels[k].children.list[0]].mom = k;
  voxels[voxels[k].children.list[1]].mom = k;
      
  delete [] pt;
    
  return;
  }
  
  
void split_tree::create_voxels(int in, double RMT[3][3], double smn, double smx, double ar, double hi[3], double lo[3], double mtol)
  {
  int code = 0, m, n; 
  double mag, tmp[3], M[3];
  Point cpt0, cpt1, p0, p1, cent;
  Point pt[6]; //will do +8 to get actual pt
  Point pnt[6]; //will do +8 to get actual pt
  Vector vec[6]; //will do +8 to get actual pt
  int nbr[6];
  double del[3];
  Vector v1;
  double begx, begy, begz, endx, endy, endz, halfx, halfy, halfz, tol, artest;

  //will loop through until facet visited all voxels, at finest level, no more refining 
  //set cornerpts
  cpt0 = voxels[in].cornerpts[0];
  cpt1 = voxels[in].cornerpts[1];
        
  //determine metric lengths
  for (m = 0; m < 3; m++)
    {
    switch (m) 
      {
      case 0:
        v1 = Vector(fabs(cpt1[0] - cpt0[0]),0.0,0.0); 
      break;
      case 1:
        v1 = Vector(0.0,fabs(cpt1[1] - cpt0[1]),0.0);
      break;
      case 2:
        v1 = Vector(0.0,0.0,fabs(cpt1[2] - cpt0[2]));
      break;
      default:
        fprintf(stderr,"\nVoxel %d has more than three dimensions.  Exiting....\n",in);
        fflush(stderr);
        //MPI_Abort(MPI_COMM_WORLD,0);
        exit(0);
      break;
      }
    //fill in M
    tmp[0]=tmp[1]=tmp[2]=0.0;
    for (n = 0; n < 3; n++)
      tmp[n]=v1[0]*RMT[0][n]+v1[1]*RMT[1][n]+v1[2]*RMT[2][n];
    mag = tmp[0]*v1[0]+tmp[1]*v1[1]+tmp[2]*v1[2];
    M[m] = sqrt(fabs(mag));
    }
    
  //now, determine if M needs changed
  //reset v1 for testing smn, smx
  v1 = Vector(fabs(cpt1[0] - cpt0[0]),fabs(cpt1[1] - cpt0[1]),fabs(cpt1[2] - cpt0[2]));
  
  //set up tol for refinement
  int mettol = (int) mtol;
  double ntol = mtol - mettol;

  if (!highquality)
  {
  //enforce desired quality constraints
  //cannot enforce gradation here, since limited to refining this voxel, but nabor takes good care
  //set up del
  for (m = 0; m < 3; m++)
    del[m] = fabs(cpt1[m] - cpt0[m]);
  
  //get centroid
  cent = cpt0 + cpt1;
  cent /= 2.0;
  
  //set up values of other pts from cornerpts
  begx = cpt0[0];
  begy = cpt0[1];
  begz = cpt0[2];
  endx = cpt1[0];
  endy = cpt1[1];
  endz = cpt1[2];
  halfx = (cpt0[0]+cpt1[0])/2.0;
  halfy = (cpt0[1]+cpt1[1])/2.0; 
  halfz = (cpt0[2]+cpt1[2])/2.0;
          
  //set up other pts of this voxel
  pt[0] = Point(halfx,halfy,begz);
  pt[1] = Point(halfx,begy,halfz);
  pt[2] = Point(endx,halfy,halfz);
  pt[3] = Point(halfx,endy,halfz);
  pt[4] = Point(begx,halfy,halfz);
  pt[5] = Point(halfx,halfy,endz);
  
  tol = 0.25*smn;
 
  for (m = 0; m < 6; m++)
    {
    vec[m] = Vector(cent,pt[m]); //get direction
    vec[m].normalize(); //normalize
    vec[m] *= tol; //get tol as length
    pnt[m] = pt[m] + Point(vec[m][0],vec[m][1],vec[m][2]);
    nbr[m] = nabor_search(pnt[m],0,smn,smx);
    }

  //check aspect ratio here...no need to worry about smn as never making smaller voxels smaller, just larger ones, and smn has thus been observed
  //cannot go smaller than smn if ar > 1, since smallest dimension already minimum smn, and larger must be greater
  //compute ar...noting that if set to one, will get isotropic
  //best thing is to refine in max direction
  //ok to use max since if one is 4.0 and the other 3.99, same result!
  artest = MAX(del[0],MAX(del[1],del[2]))/MIN(del[0],MIN(del[1],del[2]));
  //set tolerances for ar
  if (artest > (ar + ntol))
    {
    if (del[0] + (1.0e-10) > MAX(del[0],MAX(del[1],del[2])))
      {
      M[0] += 2.0;
      }
    if (del[1] + (1.0e-10) > MAX(del[0],MAX(del[1],del[2])))
      {
      M[1] += 2.0;
      }
    if (del[2] + (1.0e-10) > MAX(del[0],MAX(del[1],del[2])))
      {
      M[2] += 2.0;
      }
    }

  //VERY IMPORTANT:  This can still leave us with four neighbors on 1 given face (of opp sides...up to 3 total), 
  //and even a voxel between two refined voxels (to avoid crossbar)
  //just makes quality better by less big voxels in field of small for weird vol calc by flow solver 
  //check if any current pairs of neighbors refined, if not move on
  if (nbr[2] >= 0 && nbr[4] >= 0 && voxels[nbr[2]].split == voxels[nbr[4]].split && voxels[nbr[2]].split > 0)
    M[voxels[nbr[2]].split-1] += 2.0;
  if (nbr[1] >= 0 && nbr[3] >= 0 && voxels[nbr[1]].split == voxels[nbr[3]].split && voxels[nbr[1]].split > 0)
    M[voxels[nbr[1]].split-1] += 2.0;
  if (nbr[0] >= 0 && nbr[5] >= 0 && voxels[nbr[0]].split == voxels[nbr[5]].split && voxels[nbr[0]].split > 0)
    M[voxels[nbr[0]].split-1] += 2.0;
  }

  //if this puts us less than min spacing, don't refine (note, we might already be slightly below, but that is a better quandry than not getting there if smn not a power of 2)
  for (n = 0; n < 3; n++)
    if (v1[n] < smn-(1.0e-15))
      M[n] = 0.0;
        
  //if we are greater than max spacing, add 2.0 to each direction to assure refinement in that direction (or whichever is most dire)
  for (n = 0; n < 3; n++)
    if (v1[n] > smx+(1.0e-15))
      M[n] += 2.0;

  //now, set codes based on M (this will always favor x refinement first, as will ar calculation)
  //using tolerance covers cases of equality giving preference to the first large value and overriding the diff between 3.999 and 4.0
  //use within 1% instead of 0.0001%, although both work
  if (M[0] > mtol && (M[0] + ntol) > M[1] && (M[0] + ntol) > M[2])
    {
    code = 1;
    }
  else if (M[1] > mtol && (M[1] + ntol) > M[0] && (M[1] + ntol) > M[2])
    {
    code = 2;
    }
  else if (M[2] > mtol && (M[2] + ntol) > M[0] && (M[2] + ntol) > M[1])
    {
    code = 3;
    }
  else
    code = 0;
  
  //if code == 0, we can stop, as no need to go further down branch...simply will return 
  //else, don't refine that which has been refined, just head down tree  
  if (voxels[in].split > 0 && code > 0)
    {
    //now, see which voxels this facet is within the scope of, make recursive call
    for (m = 0; m < voxels[in].children.max; m++)
      {
      //set cornerpts
      p0 = voxels[voxels[in].children.list[m]].cornerpts[0];
      p1 = voxels[voxels[in].children.list[m]].cornerpts[1];
            
      //check if facet is at all inside this voxel...if not, continue
      if (lo[0] > p1[0] || lo[1] > p1[1] || lo[2] > p1[2] || hi[0] < p0[0] || hi[1] < p0[1] || hi[2] < p0[2])
        continue;

      create_voxels(voxels[in].children.list[m], RMT, smn, smx, ar, hi, lo, mtol);
      }
    }
  //if we have no kids and a need to refine... 
  else if (voxels[in].split == 0 && code > 0)
    {    
    refine_voxel(nvoxels, in, code);
    
    //now, see which voxels this facet is within the scope of, make recursive call
    for (m = 0; m < voxels[in].children.max; m++)
      {
      //set cornerpts
      p0 = voxels[voxels[in].children.list[m]].cornerpts[0];
      p1 = voxels[voxels[in].children.list[m]].cornerpts[1];
            
      //check if facet is at all inside this voxel...if not, continue
      if (lo[0] > p1[0] || lo[1] > p1[1] || lo[2] > p1[2] || hi[0] < p0[0] || hi[1] < p0[1] || hi[2] < p0[2])
        continue;

      create_voxels(voxels[in].children.list[m], RMT, smn, smx, ar, hi, lo, mtol);
      }
    }
    
  return;
  }

void split_tree::create_tree(geometry *geom, double &smn, double &smx, double ar, double mtol, double nspace)
  {
  int i, j, k, n, m, n0, n1, n2, flag, vt = 0;
  Vector v0, v1, v2, norm;
  double RMT[3][3], RI[3][3], Rl[3][3], lam[3][3];
  double h1, h2, h3;
  double hi[3],lo[3];
  Point p0, p1, cent;
  double begx, begy, begz, endx, endy, endz, halfx, halfy, halfz;
  double tol, delx, dely, delz, delxn, delyn, delzn, artest, mins, maxs;
  time_t tm;
  char *t_char;

  time(&tm);
  t_char = ctime(&tm);
  fprintf(out_f,"\nBeginning voxel creation at %s.\n",t_char);

  //be sure voxel_dim and nvoxels set to 0
  voxel_dim = nvoxels = 0;
  voxels = 0; //nullify voxel pointer to be sure
  
  //to set up first voxel
  vt = 1;
  
  //init voxel mem
  if (vt > voxel_dim)
    {
    my_mem(voxel_dim, nvoxels, &voxels, vt, STCHUNK);
    //set mom to -1, split/cut to 0
    for (i = 0; i < voxel_dim; i++)
      {
      voxels[i].split = 0;
      voxels[i].mom = -1;
      voxels[i].cut = 0;
      voxels[i].children.max = 0;
      for (j = 0; j < 27; j++)
        voxels[i].nodes[j] = -1;
      }
    }
  
  //now, create root cell
  voxels[0].cornerpts[0] = geom->min_point();
  voxels[0].cornerpts[1] = geom->max_point();
  
  //be sure mom set to -1
  voxels[0].mom = -1;
  //init split
  voxels[0].split = 0;
  
  //set n_voxels appropriately
  nvoxels = vt;
  
  if (!ugeom)
  {
  //loop through geom facets, create tensors, determine split
  for (i=1; i <= geom->ngb; i++)
    {
    for (j=geom->n_begin[i]; j < geom->n_end[i]; j++)
      {
      n0 = geom->g_facet[j].nodes[0];
      n1 = geom->g_facet[j].nodes[1];
      n2 = geom->g_facet[j].nodes[2];
      
      //find normal vector for any tensor formulation
      v0 = Vector(geom->g_vert[n0],geom->g_vert[n1]);
      v1 = Vector(geom->g_vert[n0],geom->g_vert[n2]);
      norm = v0 % v1;
      norm.normalize();
        
      //set up for uniform case if desired  
      if (nspace < -1.0e-06)
        {
        //find other orthogonal vector (will use v0 for second orthogonal vector)
        v2 = v0 % norm;
      
        //now, we have three orthogonal vectors, let's normalize them (will use gspace and nspace from geom for h's)
        v0.normalize();
        v2.normalize();
      
        h1 = geom->g_facet[j].space;
        h2 = geom->g_facet[j].space;
        h3 = geom->g_facet[j].n_space;
      
        //construct tensor
        RMT[0][0] = v0[0];
        RMT[1][0] = v0[1];
        RMT[2][0] = v0[2];
        RMT[0][1] = v2[0];
        RMT[1][1] = v2[1];
        RMT[2][1] = v2[2];
        RMT[0][2] = norm[0];
        RMT[1][2] = norm[1];
        RMT[2][2] = norm[2]; 
        RI[0][0] = v0[0];
        RI[0][1] = v0[1];
        RI[0][2] = v0[2];
        RI[1][0] = v2[0];
        RI[1][1] = v2[1];
        RI[1][2] = v2[2];
        RI[2][0] = norm[0];
        RI[2][1] = norm[1];
        RI[2][2] = norm[2];
        lam[0][0] = 1.0/h1/h1;
        lam[0][1] = 0.0;
        lam[0][2] = 0.0;
        lam[1][0] = 0.0;
        lam[1][1] = 1.0/h2/h2;
        lam[1][2] = 0.0;
        lam[2][0] = 0.0;
        lam[2][1] = 0.0;
        lam[2][2] = 1.0/h3/h3;

        for (n = 0; n < 3; n++)
          for (m = 0; m < 3; m++)
            Rl[m][n] = RMT[m][0]*lam[0][n]+RMT[m][1]*lam[1][n]+RMT[m][2]*lam[2][n];
        for (n = 0; n < 3; n++)
          for (m = 0; m < 3; m++)
            RMT[m][n] = Rl[m][0]*RI[0][n]+Rl[m][1]*RI[1][n]+Rl[m][2]*RI[2][n];
        }
      else
        {
        //set up specific tensor for each facet, using multiple options for normal spacing
        if (geom->layers[i] == 0)
          {
          //create centroid
          cent = (geom->g_vert[n0] + geom->g_vert[n1] + geom->g_vert[n2])/3.0;
          //now, make normal vector right length
          norm *= nspace;
          //now, make p0 (the apex of tet)
          p0 = cent + Point(norm[0],norm[1],norm[2]);
          //finally, compute tensor
          compute_Riemannian_metric(geom->g_vert[n0], geom->g_vert[n1], geom->g_vert[n2], p0, RMT);
          }
        else if (geom->layers[i] == 1)
          {
          //create centroid
          cent = (geom->g_vert[n0] + geom->g_vert[n1] + geom->g_vert[n2])/3.0;
          //now, make normal vector right length
          norm *= geom->g_facet[j].space; //using curvature
          //now, make p0 (the apex of tet)
          p0 = cent + Point(norm[0],norm[1],norm[2]);
          //finally, compute tensor
          compute_Riemannian_metric(geom->g_vert[n0], geom->g_vert[n1], geom->g_vert[n2], p0, RMT);
          }
        else if (geom->layers[i] == 2)
          {
          //create centroid
          cent = (geom->g_vert[n0] + geom->g_vert[n1] + geom->g_vert[n2])/3.0;
          //now, make normal vector right length
          norm *= geom->g_facet[j].n_space; //using intersection
          //now, make p0 (the apex of tet)
          p0 = cent + Point(norm[0],norm[1],norm[2]);
          //finally, compute tensor
          compute_Riemannian_metric(geom->g_vert[n0], geom->g_vert[n1], geom->g_vert[n2], p0, RMT);
          }
        else if (geom->layers[i] == 3)
          {
          //create centroid
          cent = (geom->g_vert[n0] + geom->g_vert[n1] + geom->g_vert[n2])/3.0;
          //now, make normal vector right length
          norm *= MIN(geom->g_facet[j].n_space,geom->g_facet[j].space); //using min of curvature/intersection
          //now, make p0 (the apex of tet)
          p0 = cent + Point(norm[0],norm[1],norm[2]);
          //finally, compute tensor
          compute_Riemannian_metric(geom->g_vert[n0], geom->g_vert[n1], geom->g_vert[n2], p0, RMT);
          }
        else if (geom->layers[i] == 4)
          {
          //create centroid
          cent = (geom->g_vert[n0] + geom->g_vert[n1] + geom->g_vert[n2])/3.0;
          //now, make normal vector right length
          norm *= MAX(geom->g_facet[j].n_space,geom->g_facet[j].space); //using max of curvature/intersection
          //now, make p0 (the apex of tet)
          p0 = cent + Point(norm[0],norm[1],norm[2]);
          //finally, compute tensor
          compute_Riemannian_metric(geom->g_vert[n0], geom->g_vert[n1], geom->g_vert[n2], p0, RMT);
          }
        else if (geom->layers[i] == 5)
          {
          //create centroid
          cent = (geom->g_vert[n0] + geom->g_vert[n1] + geom->g_vert[n2])/3.0;
          //now, make normal vector right length
          norm *= (distance(geom->g_vert[n0],geom->g_vert[n1]) + distance(geom->g_vert[n1],geom->g_vert[n2]) + distance(geom->g_vert[n2],geom->g_vert[n0]))/3.0; //average side length
          //now, make p0 (the apex of tet)
          p0 = cent + Point(norm[0],norm[1],norm[2]);
          //finally, compute tensor
          compute_Riemannian_metric(geom->g_vert[n0], geom->g_vert[n1], geom->g_vert[n2], p0, RMT);
          }
        }
          
      //set hi, lo
      for (m = 0; m < 3; m++)
        {
        hi[m] = -1.0e+20;
        lo[m] = 1.0e+20;
        }
          
      //now, set up facet extents for use in determining where to go in children
      for (n = 0; n < 3; n++)
        {
        switch (n) 
          {
          case 0:
            k = n0;
          break;
          case 1:
            k = n1;
          break;
          case 2:
            k = n2;
          break;
          default:
            fprintf(stderr,"\nFacet %d has more than three sides.  Exiting....\n",j);
            fflush(stderr);
            //MPI_Abort(MPI_COMM_WORLD,0);
            exit(0);
          break;
          }
        for (m = 0; m < 3; m++)
          {
          hi[m] = MAX(hi[m],geom->g_vert[k][m]);
          lo[m] = MIN(lo[m],geom->g_vert[k][m]);
          }
        }
      //use function by passing in root, then letting it recursively head down tree
      create_voxels(0,RMT,smn,smx,ar,hi,lo,mtol);
      }
    }
  }
    
  //set up for adapt...already much refinement due to geom 
  if (adapt == 1)
    {
    fprintf(out_f,"\nApplying adaptive spacing field.");
    fflush(out_f);
    double left[3][3], right[3][3], ev1[3], ev2[3], ev3[3], lambda[3];

    // pass adaptive tensors down tree
    int nsp = SF_Number_of_Entries();
    fprintf(out_f,"\nNumber of entries in adaptive file = %d\n",nsp);
    fflush(out_f);
    int tenth = MAX(1,nsp/10);
    int hundredth = MAX(1,nsp/100);
    int percent = 0;
    int nvrt;
    double vert[8][3];
    for (i=0; i < nsp; i++)
      {
      if ((i+1) % hundredth == 0) fprintf(out_f,".");
      if ((i+1) % tenth == 0 || (i+1) == nsp) fprintf(out_f," %d %%\n",(percent = percent + 10));
      fflush(out_f);

      SF_Retrieve_Tensor_Item(i,RMT,nvrt,vert);
      
      for (j = 0; j < 3; j++)
        {
        lo[j] = 1.0e+20;
        hi[j] = -1.0e+20;
        }

      for (j=0; j < nvrt; j++)
        {
        lo[0] = MIN(lo[0],vert[j][0]);
        lo[1] = MIN(lo[1],vert[j][1]);
        lo[2] = MIN(lo[2],vert[j][2]);
        hi[0] = MAX(hi[0],vert[j][0]);
        hi[1] = MAX(hi[1],vert[j][1]);
        hi[2] = MAX(hi[2],vert[j][2]);
        }

      //no need to do this since observing min and max spacing in create_voxels
      /*SF_Decompose_Tensor(RMT, left, right, lambda);
      h1 = sqrt(1.0/lambda[0]);
      h2 = sqrt(1.0/lambda[1]);
      h3 = sqrt(1.0/lambda[2]);
      h1 = MAX(smn,MIN(smx,h1));
      h2 = MAX(smn,MIN(smx,h2));
      h3 = MAX(smn,MIN(smx,h3));

      ev1[0] = left[0][0];
      ev1[1] = left[1][0];
      ev1[2] = left[2][0];
      ev2[0] = left[0][1];
      ev2[1] = left[1][1];
      ev2[2] = left[2][1];
      ev3[0] = left[0][2];
      ev3[1] = left[1][2];
      ev3[2] = left[2][2];

      SF_Compute_Riemannian_Metric(ev1, ev2, ev3, h1, h2, h3, RMT);*/
      
      //use function by passing in root, then letting it recursively head down tree
      create_voxels(0,RMT,smn,smx,ar,hi,lo,mtol);
      }
      
    //now, that we are done with SF, finalize
    SF_Finalize();
    }
    
  //finally, be sure max spacing is truly observed, reset min and max
  flag = 1;
  do 
    {
    mins = 1.0e+20;
    maxs = -1.0e+20; 
    flag = 0;
    for (k = 0; k < nvoxels; k++)
      {
      if (voxels[k].children.max > 0)
        continue;
        
      p0 = voxels[k].cornerpts[0];
      p1 = voxels[k].cornerpts[1];
      
      delx = fabs(p1[0] - p0[0]);
      dely = fabs(p1[1] - p0[1]);
      delz = fabs(p1[2] - p0[2]);
     
      mins = MIN(mins,MIN(delx,MIN(dely,delz)));
      maxs = MAX(maxs,MAX(delx,MAX(dely,delz)));
      
      if (delx > smx+(1.0e-12))
        {
        refine_voxel(nvoxels, k, 1);
        flag++;
        }
      else if (dely > smx+(1.0e-12))
        {
        refine_voxel(nvoxels, k, 2);
        flag++;
        }
      else if (delz > smx+(1.0e-12))
        {
        refine_voxel(nvoxels, k, 3);
        flag++;
        }

      //set up tol for refinement
      int mettol = (int) mtol;
      double ntol = mtol - mettol;

      //now, check AR
      double artest = MAX(delx,MAX(dely,delz))/MIN(delx,MIN(dely,delz));
      //set tolerances for ar
      if (artest > (ar + ntol))
        {
        if (delx + (1.0e-10) > MAX(delx,MAX(dely,delz)))
          {
          refine_voxel(nvoxels, k, 1);
          flag++;
          }
        else if (dely + (1.0e-10) > MAX(delx,MAX(dely,delz)))
          {
          refine_voxel(nvoxels, k, 2);
          flag++;
          }
        else if (delz + (1.0e-10) > MAX(delx,MAX(dely,delz)))
          {
          refine_voxel(nvoxels, k, 3);
          flag++;
          }
        }
      }
    }while (flag > 0);
  
  //now, reset min and max for neighbor routine
  if (mins < smn-(1.0e-12))
    smn = mins;
  if (maxs < smx-(1.0e-12))
    smx = maxs;
  
  time(&tm);
  t_char = ctime(&tm);
  fprintf(out_f,"\nFinished voxel creation at %s, beginning quality control.\n",t_char);

  fprintf(out_f, "\n Actual mesh spacing before quality metrics:  min = %lf, max = %lf, smn = %lf, smx = %lf\n",mins,maxs,smn,smx);
  fprintf(out_f,"\nAfter initial refinement, nvoxels = %d\n",nvoxels);
  fflush(out_f); 
        
  /*for (i = 0; i < nvoxels; i++)
    if (voxels[i].cornerpts[0][0] >= voxels[i].cornerpts[1][0] || voxels[i].cornerpts[0][1] >= voxels[i].cornerpts[1][1] || voxels[i].cornerpts[0][2] >= voxels[i].cornerpts[1][2])
      fprintf(out_f,"\n voxel %d has wrongly oriented cpts %lf %lf %lf and %lf %lf %lf\n",i,voxels[i].cornerpts[0][0],voxels[i].cornerpts[0][1],voxels[i].cornerpts[0][2],voxels[i].cornerpts[1][0],voxels[i].cornerpts[1][1],voxels[i].cornerpts[1][2]);
  fflush(out_f);*/
    
  //now, tree is created, do metric quality checks
  //look for 4-1 on anedge, then for crossbars and 4-1 on faces, then look to keep gradation to 2-1, then look to be sure ar is observed
  int total_changes = 1; //set to go through loop
  int changes1 = 0, changes2 = 0, changes3 = 0, changes4 = 0, changes5 = 0, pass = 0;
  Point *pt = new Point[18]; //will do +8 to get actual pt
  Point *pnt = new Point[18]; //will do +8 to get actual pt
  Vector *vec = new Vector[18]; //will do +8 to get actual pt
  Vector pert; //for face pert
  int fnbr[4]; 
  int nbr[18];
  List **nbrs = new List*[6];
  for (k = 0; k < 6; k++)
    nbrs[k] = new List();
  Point pt0, pt1, pt2, pt3, pt4, pt5, pt6, pt7, temp;
  flag = 0; //reset
  double del[3];
  int type[2];
  //set up index array to make bias consistent
  int inx[18];
  inx[0] = 14;
  inx[1] = 16;
  inx[2] = 15;
  inx[3] = 13;
  inx[4] = 17;
  inx[5] = 12;
  inx[6] = 0;
  inx[7] = 2;
  inx[8] = 8;
  inx[9] = 10;
  inx[10] = 1;
  inx[11] = 3;
  inx[12] = 9;
  inx[13] = 11;
  inx[14] = 4;
  inx[15] = 5;
  inx[16] = 6;
  inx[17] = 7;
  
  //set up tol for AR check
  int mettol = (int) mtol;
  double ntol = mtol - mettol;
  
  do
    {
    total_changes = 0; //reset
    changes1 = changes2 = changes3 = changes4 = changes5 = 0;
      
    for (i = 0; i < nvoxels && highquality; i++)
      {
      if (voxels[i].children.max > 0)
        continue; //only check on finest level
    
      //we need to check all directions
      //get crnerpts
      p0 = voxels[i].cornerpts[0];
      p1 = voxels[i].cornerpts[1];
      
      //set up del
      for (j = 0; j < 3; j++)
        del[j] = fabs(p1[j] - p0[j]);
  
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
      pt[12] = Point(halfx,halfy,begz);
      pt[13] = Point(halfx,begy,halfz);
      pt[14] = Point(endx,halfy,halfz);
      pt[15] = Point(halfx,endy,halfz);
      pt[16] = Point(begx,halfy,halfz);
      pt[17] = Point(halfx,halfy,endz);
      
      //set up cornerpts for use in face nabors
      pt0 = p0;
      pt1 = Point(endx,begy,begz);
      pt2 = Point(endx,endy,begz);
      pt3 = Point(begx,endy,begz);
      pt4 = Point(begx,begy,endz);
      pt5 = Point(endx,begy,endz);
      pt6 = p1;
      pt7 = Point(begx,endy,endz);
      
      //reset nbrs
      for (k = 0; k < 6; k++)
        nbrs[k]->Redimension(0);
  
      tol = 0.25*smn;
      
      flag = 0;    
      for (m = 0; m < 6 && !flag; m++)
        {
          j = inx[m];
          vec[j] = Vector(cent,pt[j]); //get direction
          vec[j].normalize(); //normalize
          vec[j] *= tol; //get tol as length
          pnt[j] = pt[j] + Point(vec[j][0],vec[j][1],vec[j][2]);
          
          //now, branch out pnt four directions
          if (j == 12)
            {
            for (k = 0; k < 4; k++)
              {
              switch (k) 
                {
                case 0:
                  pert = Vector(pt[j],pt0);
                break;
                case 1:
                  pert = Vector(pt[j],pt1);
                break;
                case 2:
                  pert = Vector(pt[j],pt2);
                break;
                case 3:
                  pert = Vector(pt[j],pt3);
                break;
                default:
                break;
                }
              pert.normalize();
              pert *= tol;
              temp = pnt[j] + Point(pert[0],pert[1],pert[2]);
              fnbr[k] = nabor_search(temp,0,smn,smx);
              }
            }
          else if (j == 13)
            {
            for (k = 0; k < 4; k++)
              {
              switch (k) 
                {
                case 0:
                  pert = Vector(pt[j],pt0);
                break;
                case 1:
                  pert = Vector(pt[j],pt1);
                break;
                case 2:
                  pert = Vector(pt[j],pt5);
                break;
                case 3:
                  pert = Vector(pt[j],pt4);
                break;
                default:
                break;
                }
              pert.normalize();
              pert *= tol;
              temp = pnt[j] + Point(pert[0],pert[1],pert[2]);
              fnbr[k] = nabor_search(temp,0,smn,smx);
              }
            }
          else if (j == 14)
            {
            for (k = 0; k < 4; k++)
              {
              switch (k) 
                {
                case 0:
                  pert = Vector(pt[j],pt1);
                break;
                case 1:
                  pert = Vector(pt[j],pt2);
                break;
                case 2:
                  pert = Vector(pt[j],pt6);
                break;
                case 3:
                  pert = Vector(pt[j],pt5);
                break;
                default:
                break;
                }
              pert.normalize();
              pert *= tol;
              temp = pnt[j] + Point(pert[0],pert[1],pert[2]);
              fnbr[k] = nabor_search(temp,0,smn,smx);
              }
            }
          else if (j == 15)
            {
            for (k = 0; k < 4; k++)
              {
              switch (k) 
                {
                case 0:
                  pert = Vector(pt[j],pt2);
                break;
                case 1:
                  pert = Vector(pt[j],pt3);
                break;
                case 2:
                  pert = Vector(pt[j],pt7);
                break;
                case 3:
                  pert = Vector(pt[j],pt6);
                break;
                default:
                break;
                }
              pert.normalize();
              pert *= tol;
              temp = pnt[j] + Point(pert[0],pert[1],pert[2]);
              fnbr[k] = nabor_search(temp,0,smn,smx);
              }
            }
          else if (j == 16)
            {
            for (k = 0; k < 4; k++)
              {
              switch (k) 
                {
                case 0:
                  pert = Vector(pt[j],pt3);
                break;
                case 1:
                  pert = Vector(pt[j],pt0);
                break;
                case 2:
                  pert = Vector(pt[j],pt4);
                break;
                case 3:
                  pert = Vector(pt[j],pt7);
                break;
                default:
                break;
                }
              pert.normalize();
              pert *= tol;
              temp = pnt[j] + Point(pert[0],pert[1],pert[2]);
              fnbr[k] = nabor_search(temp,0,smn,smx);
              }
            }
          else if (j == 17)
            {
            for (k = 0; k < 4; k++)
              {
              switch (k) 
                {
                case 0:
                  pert = Vector(pt[j],pt6);
                break;
                case 1:
                  pert = Vector(pt[j],pt7);
                break;
                case 2:
                  pert = Vector(pt[j],pt4);
                break;
                case 3:
                  pert = Vector(pt[j],pt5);
                break;
                default:
                break;
                }
              pert.normalize();
              pert *= tol;
              temp = pnt[j] + Point(pert[0],pert[1],pert[2]);
              fnbr[k] = nabor_search(temp,0,smn,smx);
              }
            }
          //now, set up proper amount of nbrs to test
          for (k = 0; k < 4; k++)
            nbrs[j-12]->Check_List(fnbr[k]);
        
        if (nbrs[j-12]->max == 1 && nbrs[j-12]->list[0] == -1)
          continue; //nothing to worry about
        
        //set up for nabor deltas        
        int deln[3];                
        
        //looking at z-faces
        if (j == 12 || j == 17)
          {
          for (k = 0; k < nbrs[j-12]->max && !flag; k++)
            {
            nbr[j] = nbrs[j-12]->list[k];
            
            //bump out if non-existant
            if (nbr[j] < 0)
              continue;
            
            //check for errant nabor
            if (voxels[nbr[j]].children.max > 0)
              {
              fprintf(out_f,"\nYou have found an errant face %d neighbor = %d of voxel %d with kids!  Exiting....\n",j-12,nbr[j],i);
              fflush(out_f);
              voxels[nbr[j]].print(out_f);
              fflush(out_f);
              voxels[i].print(out_f);
              fflush(out_f);
              exit(0);
              }    
            
            //find neighbor delta
            deln[2] = fabs(voxels[nbr[j]].cornerpts[1][2] - voxels[nbr[j]].cornerpts[0][2]);
            
            if (deln[2] > 3.0*del[2])
              {
              //do gradation changes
              refine_voxel(nvoxels, nbr[j], 3);
              changes1++;
              //flag = 1;
              }                        
            }
          }
        //looking at y-faces
        if (j == 13 || j == 15)
          {          
          for (k = 0; k < nbrs[j-12]->max && !flag; k++)
            {
            nbr[j] = nbrs[j-12]->list[k];
            
            //bump out if non-existant
            if (nbr[j] < 0)
              continue;
            
            //check for errant nabor
            if (voxels[nbr[j]].children.max > 0)
              {
              fprintf(out_f,"\nYou have found an errant face %d neighbor = %d of voxel %d with kids!  Exiting....\n",j-12,nbr[j],i);
              fflush(out_f);
              voxels[nbr[j]].print(out_f);
              fflush(out_f);
              voxels[i].print(out_f);
              fflush(out_f);
              exit(0);
              }               
            
            //find neighbor delta
            deln[1] = fabs(voxels[nbr[j]].cornerpts[1][1] - voxels[nbr[j]].cornerpts[0][1]);
                                                
            if (deln[1] > 3.0*del[1])
              {
              //do gradation changes
              refine_voxel(nvoxels, nbr[j], 2);
              changes1++;
              //flag = 1;
              }
            }
          }
        //looking at x-faces
        if (j == 14 || j == 16)
          {          
          for (k = 0; k < nbrs[j-12]->max && !flag; k++)
            {
            nbr[j] = nbrs[j-12]->list[k];
            
            //bump out if non-existant
            if (nbr[j] < 0)
              continue;
            
            //check for errant nabor
            if (voxels[nbr[j]].children.max > 0)
              {
              fprintf(out_f,"\nYou have found an errant face %d neighbor = %d of voxel %d with kids!  Exiting....\n",j-12,nbr[j],i);
              fflush(out_f);
              voxels[nbr[j]].print(out_f);
              fflush(out_f);
              voxels[i].print(out_f);
              fflush(out_f);
              exit(0);
              }
              
            //find neighbor delta
            deln[0] = fabs(voxels[nbr[j]].cornerpts[1][0] - voxels[nbr[j]].cornerpts[0][0]);
                        
            if (deln[0] > 3.0*del[0])
              {
              //do gradation changes
              refine_voxel(nvoxels, nbr[j], 1);
              changes1++;
              //flag = 1;
              }
            }
          }
        }
        
      if (flag)
        {
        //just to be sure smn still meaningful
        for (k = 0; k < 3; k++)
          if (del[k]/2.0 < smn)
            smn = del[k]/2.0;
            
        continue; //don't do nabor gradation if already refined
        }
      
      //VERY IMPORTANT:  This can still leave us with four neighbors on 1 given face (of opp sides...up to 3 total), 
      //and even a voxel between two refined voxels (to avoid crossbar)
      //just makes quality better by less big voxels in field of small for weird vol calc by flow solver 
      //check if any current pairs of neighbors refined, if not move on
      if (!flag && nbrs[2]->max > 1 && nbrs[4]->max > 1)
        flag = voxel_grad(nbrs[2], nbrs[4], 1, 2, del, i, changes2, smn, smx, pnt[14], pnt[16]);
      else if (!flag && nbrs[1]->max > 1 && nbrs[3]->max > 1)
        flag = voxel_grad(nbrs[1], nbrs[3], 2, 0, del, i, changes2, smn, smx, pnt[13], pnt[15]);
      else if (!flag && nbrs[0]->max > 1 && nbrs[5]->max > 1)
        flag = voxel_grad(nbrs[0], nbrs[5], 0, 1, del, i, changes2, smn, smx, pnt[12], pnt[17]);
      
      if (flag)
        {
        //just to be sure smn still meaningful
        for (k = 0; k < 3; k++)
          if (del[k]/2.0 < smn)
            smn = del[k]/2.0;
            
        continue; //don't do AR if already refined
        }

      //check aspect ratio here...no need to worry about smn as never making smaller voxels smaller, just larger ones, and smn has thus been observed
      //cannot go smaller than smn if ar > 1, since smallest dimension already minimum smn, and larger must be greater
      //compute ar...noting that if set to one, will get isotropic
      //best thing is to refine in max direction
      //ok to use max since if one is 4.0 and the other 3.99, same result!
      artest = MAX(del[0],MAX(del[1],del[2]))/MIN(del[0],MIN(del[1],del[2]));
      //set tolerances for ar
      if (artest > (ar + ntol))
        {
        if (del[0] + (1.0e-10) > MAX(del[0],MAX(del[1],del[2])))
          {
          refine_voxel(nvoxels,i, 1);
          changes3++;
          }
        else if (del[1] + (1.0e-10) > MAX(del[0],MAX(del[1],del[2])))
          {
          refine_voxel(nvoxels,i, 2);
          changes3++;
          }
        else if (del[2] + (1.0e-10) > MAX(del[0],MAX(del[1],del[2])))
          {
          refine_voxel(nvoxels,i, 3);
          changes3++;
          }
        }
        
      //just to be sure smn still meaningful
      for (k = 0; k < 3; k++)
        if (del[k]/2.0 < smn)
          smn = del[k]/2.0;
      }
      
    for (i = 0; i < nvoxels; i++)
      {
      if (voxels[i].children.max > 0)
        continue; //only check on finest level
    
      //we need to check all directions
      //get crnerpts
      p0 = voxels[i].cornerpts[0];
      p1 = voxels[i].cornerpts[1];
      
      //set up del
      for (j = 0; j < 3; j++)
        del[j] = fabs(p1[j] - p0[j]);
  
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
      pt[0] = Point(halfx,begy,begz);
      pt[1] = Point(endx,halfy,begz);
      pt[2] = Point(halfx,endy,begz);
      pt[3] = Point(begx,halfy,begz);
      pt[4] = Point(begx,begy,halfz);
      pt[5] = Point(endx,begy,halfz);
      pt[6] = Point(endx,endy,halfz);
      pt[7] = Point(begx,endy,halfz);
      pt[8] = Point(halfx,begy,endz);
      pt[9] = Point(endx,halfy,endz);
      pt[10] = Point(halfx,endy,endz);
      pt[11] = Point(begx,halfy,endz);
      pt[12] = Point(halfx,halfy,begz);
      pt[13] = Point(halfx,begy,halfz);
      pt[14] = Point(endx,halfy,halfz);
      pt[15] = Point(halfx,endy,halfz);
      pt[16] = Point(begx,halfy,halfz);
      pt[17] = Point(halfx,halfy,endz);
      
      //set up cornerpts for use in face nabors
      pt0 = p0;
      pt1 = Point(endx,begy,begz);
      pt2 = Point(endx,endy,begz);
      pt3 = Point(begx,endy,begz);
      pt4 = Point(begx,begy,endz);
      pt5 = Point(endx,begy,endz);
      pt6 = p1;
      pt7 = Point(begx,endy,endz);
  
      tol = 0.25*smn;
      
      //reset nbrs
      for (k = 0; k < 6; k++)
        nbrs[k]->Redimension(0);
      
      //since we should NEVER have less than 0.5 smn as a voxel size, 
      //the pert of .25smn each direction for face perts yields 1/(2sqrt2) 
      //smn hypotenuse, which should be well within any voxel tested!
      //do CB first, since any refinement here might mitigate the need for 4-1 elsewhere
      flag = 0;    
      for (m = 0; m < 6 && !flag; m++)
        {
        j = inx[m];
        //non-face
        if (j < 12)
          {    
          vec[j] = Vector(cent,pt[j]); //get direction
          vec[j].normalize(); //normalize
          vec[j] *= tol; //get tol as length
          pnt[j] = pt[j] + Point(vec[j][0],vec[j][1],vec[j][2]);
          nbr[j] = nabor_search(pnt[j],0,smn,smx);
          //if (nbr[j] >= 0)
            //voxels[nbr[j]].print(out_f);
          }
        else
          {
          vec[j] = Vector(cent,pt[j]); //get direction
          vec[j].normalize(); //normalize
          vec[j] *= tol; //get tol as length
          pnt[j] = pt[j] + Point(vec[j][0],vec[j][1],vec[j][2]);
          
          //now, branch out pnt four directions
          if (j == 12)
            {
            for (k = 0; k < 4; k++)
              {
              switch (k) 
                {
                case 0:
                  pert = Vector(pt[j],pt0);
                break;
                case 1:
                  pert = Vector(pt[j],pt1);
                break;
                case 2:
                  pert = Vector(pt[j],pt2);
                break;
                case 3:
                  pert = Vector(pt[j],pt3);
                break;
                default:
                break;
                }
              pert.normalize();
              pert *= tol;
              temp = pnt[j] + Point(pert[0],pert[1],pert[2]);
              fnbr[k] = nabor_search(temp,0,smn,smx);
              }
            }
          else if (j == 13)
            {
            for (k = 0; k < 4; k++)
              {
              switch (k) 
                {
                case 0:
                  pert = Vector(pt[j],pt0);
                break;
                case 1:
                  pert = Vector(pt[j],pt1);
                break;
                case 2:
                  pert = Vector(pt[j],pt5);
                break;
                case 3:
                  pert = Vector(pt[j],pt4);
                break;
                default:
                break;
                }
              pert.normalize();
              pert *= tol;
              temp = pnt[j] + Point(pert[0],pert[1],pert[2]);
              fnbr[k] = nabor_search(temp,0,smn,smx);
              }
            }
          else if (j == 14)
            {
            for (k = 0; k < 4; k++)
              {
              switch (k) 
                {
                case 0:
                  pert = Vector(pt[j],pt1);
                break;
                case 1:
                  pert = Vector(pt[j],pt2);
                break;
                case 2:
                  pert = Vector(pt[j],pt6);
                break;
                case 3:
                  pert = Vector(pt[j],pt5);
                break;
                default:
                break;
                }
              pert.normalize();
              pert *= tol;
              temp = pnt[j] + Point(pert[0],pert[1],pert[2]);
              fnbr[k] = nabor_search(temp,0,smn,smx);
              }
            }
          else if (j == 15)
            {
            for (k = 0; k < 4; k++)
              {
              switch (k) 
                {
                case 0:
                  pert = Vector(pt[j],pt2);
                break;
                case 1:
                  pert = Vector(pt[j],pt3);
                break;
                case 2:
                  pert = Vector(pt[j],pt7);
                break;
                case 3:
                  pert = Vector(pt[j],pt6);
                break;
                default:
                break;
                }
              pert.normalize();
              pert *= tol;
              temp = pnt[j] + Point(pert[0],pert[1],pert[2]);
              fnbr[k] = nabor_search(temp,0,smn,smx);
              }
            }
          else if (j == 16)
            {
            for (k = 0; k < 4; k++)
              {
              switch (k) 
                {
                case 0:
                  pert = Vector(pt[j],pt3);
                break;
                case 1:
                  pert = Vector(pt[j],pt0);
                break;
                case 2:
                  pert = Vector(pt[j],pt4);
                break;
                case 3:
                  pert = Vector(pt[j],pt7);
                break;
                default:
                break;
                }
              pert.normalize();
              pert *= tol;
              temp = pnt[j] + Point(pert[0],pert[1],pert[2]);
              fnbr[k] = nabor_search(temp,0,smn,smx);
              }
            }
          else if (j == 17)
            {
            for (k = 0; k < 4; k++)
              {
              switch (k) 
                {
                case 0:
                  pert = Vector(pt[j],pt6);
                break;
                case 1:
                  pert = Vector(pt[j],pt7);
                break;
                case 2:
                  pert = Vector(pt[j],pt4);
                break;
                case 3:
                  pert = Vector(pt[j],pt5);
                break;
                default:
                break;
                }
              pert.normalize();
              pert *= tol;
              temp = pnt[j] + Point(pert[0],pert[1],pert[2]);
              fnbr[k] = nabor_search(temp,0,smn,smx);
              }
            }
          //now, set up proper amount of nbrs to test
          for (k = 0; k < 4; k++)
            nbrs[j-12]->Check_List(fnbr[k]);
          }
        
        if (nbrs[j-12]->max == 1 && nbrs[j-12]->list[0] == -1)
          continue; //nothing to worry about
        
        //looking at z-faces
        if (j == 12 || j == 17)
          {
          //set type
          type[0] = 1;
          type[1] = 2;
            
          for (k = 0; k < nbrs[j-12]->max && !flag; k++)
            {
            nbr[j] = nbrs[j-12]->list[k];
            
            //bump out if non-existant
            if (nbr[j] < 0)
              continue;
            
            //check for errant nabor
            if (voxels[nbr[j]].children.max > 0)
              {
              fprintf(out_f,"\nYou have found an errant face %d neighbor = %d of voxel %d with kids!  Exiting....\n",j-12,nbr[j],i);
              fflush(out_f);
              voxels[nbr[j]].print(out_f);
              fflush(out_f);
              voxels[i].print(out_f);
              fflush(out_f);
              exit(0);
              }    
                        
            flag = eval_qualityCB(i, nbr[j], del, type, changes4);
            }
          }
        //looking at y-faces
        if (j == 13 || j == 15)
          {
          //set type
          type[0] = 3;
          type[1] = 1;
          
          for (k = 0; k < nbrs[j-12]->max && !flag; k++)
            {
            nbr[j] = nbrs[j-12]->list[k];
            
            //bump out if non-existant
            if (nbr[j] < 0)
              continue;
            
            //check for errant nabor
            if (voxels[nbr[j]].children.max > 0)
              {
              fprintf(out_f,"\nYou have found an errant face %d neighbor = %d of voxel %d with kids!  Exiting....\n",j-12,nbr[j],i);
              fflush(out_f);
              voxels[nbr[j]].print(out_f);
              fflush(out_f);
              voxels[i].print(out_f);
              fflush(out_f);
              exit(0);
              }               
                        
            flag = eval_qualityCB(i, nbr[j], del, type, changes4);
            }
          }
        //looking at x-faces
        if (j == 14 || j == 16)
          { 
          //set type
          type[0] = 2;
          type[1] = 3;         
          
          for (k = 0; k < nbrs[j-12]->max && !flag; k++)
            {
            nbr[j] = nbrs[j-12]->list[k];
            
            //bump out if non-existant
            if (nbr[j] < 0)
              continue;
            
            //check for errant nabor
            if (voxels[nbr[j]].children.max > 0)
              {
              fprintf(out_f,"\nYou have found an errant face %d neighbor = %d of voxel %d with kids!  Exiting....\n",j-12,nbr[j],i);
              fflush(out_f);
              voxels[nbr[j]].print(out_f);
              fflush(out_f);
              voxels[i].print(out_f);
              fflush(out_f);
              exit(0);
              }  
                        
            flag = eval_qualityCB(i, nbr[j], del, type, changes4);
            }
          }
        }      
        
      //just to be sure smn still meaningful
      for (k = 0; k < 3; k++)
        if (del[k]/2.0 < smn)
          smn = del[k]/2.0;
      }
      
    for (i = 0; i < nvoxels; i++)
      {
      if (voxels[i].children.max > 0)
        continue; //only check on finest level
    
      //we need to check all directions
      //get crnerpts
      p0 = voxels[i].cornerpts[0];
      p1 = voxels[i].cornerpts[1];
      
      //set up del
      for (j = 0; j < 3; j++)
        del[j] = fabs(p1[j] - p0[j]);
  
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
      pt[0] = Point(halfx,begy,begz);
      pt[1] = Point(endx,halfy,begz);
      pt[2] = Point(halfx,endy,begz);
      pt[3] = Point(begx,halfy,begz);
      pt[4] = Point(begx,begy,halfz);
      pt[5] = Point(endx,begy,halfz);
      pt[6] = Point(endx,endy,halfz);
      pt[7] = Point(begx,endy,halfz);
      pt[8] = Point(halfx,begy,endz);
      pt[9] = Point(endx,halfy,endz);
      pt[10] = Point(halfx,endy,endz);
      pt[11] = Point(begx,halfy,endz);
      pt[12] = Point(halfx,halfy,begz);
      pt[13] = Point(halfx,begy,halfz);
      pt[14] = Point(endx,halfy,halfz);
      pt[15] = Point(halfx,endy,halfz);
      pt[16] = Point(begx,halfy,halfz);
      pt[17] = Point(halfx,halfy,endz);
      
      //set up cornerpts for use in face nabors
      pt0 = p0;
      pt1 = Point(endx,begy,begz);
      pt2 = Point(endx,endy,begz);
      pt3 = Point(begx,endy,begz);
      pt4 = Point(begx,begy,endz);
      pt5 = Point(endx,begy,endz);
      pt6 = p1;
      pt7 = Point(begx,endy,endz);
  
      tol = 0.25*smn;
      
      //reset nbrs
      for (k = 0; k < 6; k++)
        nbrs[k]->Redimension(0);
      
      //don't bump out after 4-1 to make better symmetry       
      //since we should NEVER have less than 0.5 smn as a voxel size, 
      //the pert of .25smn each direction for face perts yields 1/(2sqrt2) 
      //smn hypotenuse, which should be well within any voxel tested!
      //use cut as holder for 4-1 to avoid over-refinement due to 4-1 caddy corner
      for (m = 0; m < 18; m++)    
        {
        j = inx[m];
        //non-face
        if (j < 12)
          {    
          vec[j] = Vector(cent,pt[j]); //get direction
          vec[j].normalize(); //normalize
          vec[j] *= tol; //get tol as length
          pnt[j] = pt[j] + Point(vec[j][0],vec[j][1],vec[j][2]);
          nbr[j] = nabor_search(pnt[j],0,smn,smx);
          //if (nbr[j] >= 0)
            //voxels[nbr[j]].print(out_f);
          }
        else
          {
          vec[j] = Vector(cent,pt[j]); //get direction
          vec[j].normalize(); //normalize
          vec[j] *= tol; //get tol as length
          pnt[j] = pt[j] + Point(vec[j][0],vec[j][1],vec[j][2]);
          
          //now, branch out pnt four directions
          if (j == 12)
            {
            for (k = 0; k < 4; k++)
              {
              switch (k) 
                {
                case 0:
                  pert = Vector(pt[j],pt0);
                break;
                case 1:
                  pert = Vector(pt[j],pt1);
                break;
                case 2:
                  pert = Vector(pt[j],pt2);
                break;
                case 3:
                  pert = Vector(pt[j],pt3);
                break;
                default:
                break;
                }
              pert.normalize();
              pert *= tol;
              temp = pnt[j] + Point(pert[0],pert[1],pert[2]);
              fnbr[k] = nabor_search(temp,0,smn,smx);
              }
            }
          else if (j == 13)
            {
            for (k = 0; k < 4; k++)
              {
              switch (k) 
                {
                case 0:
                  pert = Vector(pt[j],pt0);
                break;
                case 1:
                  pert = Vector(pt[j],pt1);
                break;
                case 2:
                  pert = Vector(pt[j],pt5);
                break;
                case 3:
                  pert = Vector(pt[j],pt4);
                break;
                default:
                break;
                }
              pert.normalize();
              pert *= tol;
              temp = pnt[j] + Point(pert[0],pert[1],pert[2]);
              fnbr[k] = nabor_search(temp,0,smn,smx);
              }
            }
          else if (j == 14)
            {
            for (k = 0; k < 4; k++)
              {
              switch (k) 
                {
                case 0:
                  pert = Vector(pt[j],pt1);
                break;
                case 1:
                  pert = Vector(pt[j],pt2);
                break;
                case 2:
                  pert = Vector(pt[j],pt6);
                break;
                case 3:
                  pert = Vector(pt[j],pt5);
                break;
                default:
                break;
                }
              pert.normalize();
              pert *= tol;
              temp = pnt[j] + Point(pert[0],pert[1],pert[2]);
              fnbr[k] = nabor_search(temp,0,smn,smx);
              }
            }
          else if (j == 15)
            {
            for (k = 0; k < 4; k++)
              {
              switch (k) 
                {
                case 0:
                  pert = Vector(pt[j],pt2);
                break;
                case 1:
                  pert = Vector(pt[j],pt3);
                break;
                case 2:
                  pert = Vector(pt[j],pt7);
                break;
                case 3:
                  pert = Vector(pt[j],pt6);
                break;
                default:
                break;
                }
              pert.normalize();
              pert *= tol;
              temp = pnt[j] + Point(pert[0],pert[1],pert[2]);
              fnbr[k] = nabor_search(temp,0,smn,smx);
              }
            }
          else if (j == 16)
            {
            for (k = 0; k < 4; k++)
              {
              switch (k) 
                {
                case 0:
                  pert = Vector(pt[j],pt3);
                break;
                case 1:
                  pert = Vector(pt[j],pt0);
                break;
                case 2:
                  pert = Vector(pt[j],pt4);
                break;
                case 3:
                  pert = Vector(pt[j],pt7);
                break;
                default:
                break;
                }
              pert.normalize();
              pert *= tol;
              temp = pnt[j] + Point(pert[0],pert[1],pert[2]);
              fnbr[k] = nabor_search(temp,0,smn,smx);
              }
            }
          else if (j == 17)
            {
            for (k = 0; k < 4; k++)
              {
              switch (k) 
                {
                case 0:
                  pert = Vector(pt[j],pt6);
                break;
                case 1:
                  pert = Vector(pt[j],pt7);
                break;
                case 2:
                  pert = Vector(pt[j],pt4);
                break;
                case 3:
                  pert = Vector(pt[j],pt5);
                break;
                default:
                break;
                }
              pert.normalize();
              pert *= tol;
              temp = pnt[j] + Point(pert[0],pert[1],pert[2]);
              fnbr[k] = nabor_search(temp,0,smn,smx);
              }
            }
          //now, set up proper amount of nbrs to test
          for (k = 0; k < 4; k++)
            nbrs[j-12]->Check_List(fnbr[k]);
          }
        
        if ((j < 12 && nbr[j] < 0) || (j > 11 && nbrs[j-12]->max == 1 && nbrs[j-12]->list[0] == -1))
          continue; //nothing to worry about
        
        //calc deltas of nabor and kick out if illegit neighbor for single case
        if (j < 12)
          {
          //if all good, calc deltas
          delxn = fabs(voxels[nbr[j]].cornerpts[1][0] - voxels[nbr[j]].cornerpts[0][0]);
          delyn = fabs(voxels[nbr[j]].cornerpts[1][1] - voxels[nbr[j]].cornerpts[0][1]);
          delzn = fabs(voxels[nbr[j]].cornerpts[1][2] - voxels[nbr[j]].cornerpts[0][2]);
          }
                
        //looking at x-edges
        if (j == 0 || j == 2 || j == 8 || j == 10)
          {
          //since will always be a power of 2, go for 3 since that knocks out tolerance issues
          if (delxn > 3.0*del[0] && voxels[nbr[j]].children.max == 0)
            { 
            if (voxels[nbr[j]].cut == 0 || voxels[nbr[j]].cut > 1)
              voxels[nbr[j]].cut = 1; 
            }
          else if (delxn > 3.0*del[0] && voxels[nbr[j]].children.max > 0)
            {
            fprintf(out_f,"\nYou have found an errant 4-1 %d neighbor = %d of voxel %d with kids!  Exiting....\n",j-12,nbr[j],i);
            fflush(out_f);
            voxels[nbr[j]].print(out_f);
            fflush(out_f);
            voxels[i].print(out_f);
            fflush(out_f);
            exit(0);
            }
          }
        //looking at y-edges
        if (j == 1 || j == 3 || j == 9 || j == 11)
          {
          //since will always be a power of 2, go for 3 since that knocks out tolerance issues
          if (delyn > 3.0*del[1] && voxels[nbr[j]].children.max == 0)
            { 
            if (voxels[nbr[j]].cut == 0 || voxels[nbr[j]].cut > 2)
              voxels[nbr[j]].cut = 2;
            }
          else if (delyn > 3.0*del[1] && voxels[nbr[j]].children.max > 0)
            {
            fprintf(out_f,"\nYou have found an errant 4-1 %d neighbor = %d of voxel %d with kids!  Exiting....\n",j-12,nbr[j],i);
            fflush(out_f);
            voxels[nbr[j]].print(out_f);
            fflush(out_f);
            voxels[i].print(out_f);
            fflush(out_f);
            exit(0);
            }
          }
        //looking at z-edges
        if (j == 4 || j == 5 || j == 6 || j == 7)
          {          
          //since will always be a power of 2, go for 3 since that knocks out tolerance issues
          if (delzn > 3.0*del[2] && voxels[nbr[j]].children.max == 0)
            { 
            if (voxels[nbr[j]].cut == 0 || voxels[nbr[j]].cut > 3)
              voxels[nbr[j]].cut = 3;
            }
          else if (delzn > 3.0*del[2] && voxels[nbr[j]].children.max > 0)
            {
            fprintf(out_f,"\nYou have found an errant 4-1 %d neighbor = %d of voxel %d with kids!  Exiting....\n",j-12,nbr[j],i);
            fflush(out_f);
            voxels[nbr[j]].print(out_f);
            fflush(out_f);
            voxels[i].print(out_f);
            fflush(out_f);
            exit(0);
            }
          }
        
        //looking at z-faces
        if (j == 12 || j == 17)
          {
          //set type
          type[0] = 1;
          type[1] = 2;
            
          for (k = 0; k < nbrs[j-12]->max && !flag; k++)
            {
            nbr[j] = nbrs[j-12]->list[k];
            
            //bump out if non-existant
            if (nbr[j] < 0)
              continue;
            
            //check for errant nabor
            if (voxels[nbr[j]].children.max > 0)
              {
              fprintf(out_f,"\nYou have found an errant face %d neighbor = %d of voxel %d with kids!  Exiting....\n",j-12,nbr[j],i);
              fflush(out_f);
              voxels[nbr[j]].print(out_f);
              fflush(out_f);
              voxels[i].print(out_f);
              fflush(out_f);
              exit(0);
              }                            
            eval_quality41(nbr[j], del, type);
            }
          }
        //looking at y-faces
        if (j == 13 || j == 15)
          {
          //set type
          type[0] = 3;
          type[1] = 1;
          
          for (k = 0; k < nbrs[j-12]->max && !flag; k++)
            {
            nbr[j] = nbrs[j-12]->list[k];
            
            //bump out if non-existant
            if (nbr[j] < 0)
              continue;
            
            //check for errant nabor
            if (voxels[nbr[j]].children.max > 0)
              {
              fprintf(out_f,"\nYou have found an errant face %d neighbor = %d of voxel %d with kids!  Exiting....\n",j-12,nbr[j],i);
              fflush(out_f);
              voxels[nbr[j]].print(out_f);
              fflush(out_f);
              voxels[i].print(out_f);
              fflush(out_f);
              exit(0);
              }               
            
            eval_quality41(nbr[j], del, type);
            }
          }
        //looking at x-faces
        if (j == 14 || j == 16)
          { 
          //set type
          type[0] = 2;
          type[1] = 3;         
          
          for (k = 0; k < nbrs[j-12]->max && !flag; k++)
            {
            nbr[j] = nbrs[j-12]->list[k];
            
            //bump out if non-existant
            if (nbr[j] < 0)
              continue;
            
            //check for errant nabor
            if (voxels[nbr[j]].children.max > 0)
              {
              fprintf(out_f,"\nYou have found an errant face %d neighbor = %d of voxel %d with kids!  Exiting....\n",j-12,nbr[j],i);
              fflush(out_f);
              voxels[nbr[j]].print(out_f);
              fflush(out_f);
              voxels[i].print(out_f);
              fflush(out_f);
              exit(0);
              }  

            eval_quality41(nbr[j], del, type);
            }
          }
        }      
      }
      
    //now, actually refine voxels and reset cut
    //this way, we only loop through twice instead of three times...or four since this was combined with crossbars
    for (i = 0; i < nvoxels; i++)
      {
      if (voxels[i].cut > 0)
        {
        if (voxels[i].children.max == 0)
          {
          refine_voxel(nvoxels, i, voxels[i].cut);
          changes5++;
          voxels[i].cut = 0;
          }
        else
          voxels[i].cut = 0;
        }
      //set up del just to be sure smn still meaningful
      for (k = 0; k < 3; k++)
        {
        del[k] = fabs(voxels[i].cornerpts[1][k] - voxels[i].cornerpts[0][k]);
        if (del[k]/2.0 < smn)
          smn = del[k]/2.0;
        }
      }

      //finally, be sure max spacing is truly observed, reset min and max
      flag = 1;
      do 
        {
        mins = 1.0e+20;
        maxs = -1.0e+20; 
        flag = 0;
        for (k = 0; k < nvoxels; k++)
          {
          if (voxels[k].children.max > 0)
            continue;
        
          p0 = voxels[k].cornerpts[0];
          p1 = voxels[k].cornerpts[1];
      
          delx = fabs(p1[0] - p0[0]);
          dely = fabs(p1[1] - p0[1]);
          delz = fabs(p1[2] - p0[2]);
     
          mins = MIN(mins,MIN(delx,MIN(dely,delz)));
          maxs = MAX(maxs,MAX(delx,MAX(dely,delz)));
      
         if (delx > smx+(1.0e-12))
          {
          refine_voxel(nvoxels, k, 1);
          flag++;
          changes3++;
          }
        else if (dely > smx+(1.0e-12))
          {
          refine_voxel(nvoxels, k, 2);
          flag++;
          changes3++;
          }
        else if (delz > smx+(1.0e-12))
          {
          refine_voxel(nvoxels, k, 3);
          flag++;
          changes3++;
          }

        //set up tol for refinement
        int mettol = (int) mtol;
        double ntol = mtol - mettol;

        //now, check AR
        double artest = MAX(delx,MAX(dely,delz))/MIN(delx,MIN(dely,delz));
        //set tolerances for ar
        if (artest > (ar + ntol))
          {
          if (delx + (1.0e-10) > MAX(delx,MAX(dely,delz)))
            {
            refine_voxel(nvoxels, k, 1);
            flag++;
            changes3++;
            }
          else if (dely + (1.0e-10) > MAX(delx,MAX(dely,delz)))
            {
            refine_voxel(nvoxels, k, 2);
            flag++;
            changes3++;
            }
          else if (delz + (1.0e-10) > MAX(delx,MAX(dely,delz)))
            {
            refine_voxel(nvoxels, k, 3);
            flag++;
            changes3++;
            }
          }
        }
      }while (flag > 0);
  
    //now, reset min and max for neighbor routine
    if (mins < smn-(1.0e-12))
      smn = mins;
    if (maxs < smx-(1.0e-12))
      smx = maxs;
   
    //set up to avoid infinite loop
    pass++;
    
    //tally and print all changes
    total_changes = changes1+changes2+changes3+changes4+changes5;
    
    fprintf(out_f,"\nMesh quality metrics loop %d, changes (%d, %d, %d, %d, %d) = %d\n",pass,changes1,changes2,changes3,changes4,changes5,total_changes);
    fflush(out_f);
    
    }while(total_changes > 0 && pass < 1000);
  
  //delete current ptrs, regardless
  for (k = 0; k < 6; k++)
    delete nbrs[k];
  delete [] nbrs; 
  
  /*//check spacing again!
  mins = 1.0e+20;
  maxs = -1.0e+20;  
  for (k = 0; k < nvoxels; k++)
    {
    if (voxels[k].children.max > 0)
      continue;
    p0 = voxels[k].cornerpts[0];
    p1 = voxels[k].cornerpts[1];
      
    delx = fabs(p1[0] - p0[0]);
    dely = fabs(p1[1] - p0[1]);
    delz = fabs(p1[2] - p0[2]);
    
    mins = MIN(mins,MIN(delx,MIN(dely,delz)));
    maxs = MAX(maxs,MAX(delx,MAX(dely,delz)));
    }
  
  //reset for cutting
  if (mins < smn-(1.0e-12))
    smn = mins;
  if (maxs < smx-(1.0e-12))
    smx = maxs;*/

  time(&tm);
  t_char = ctime(&tm);
  fprintf(out_f,"\nFinished quality control at %s.\n",t_char);
  
  fprintf(out_f, "\n Actual mesh spacing after quality metrics:  min = %lf, max = %lf, smn = %lf, smx = %lf\n",mins,maxs,smn,smx);
  fprintf(out_f,"\nAfter mesh quality refinement, nvoxels = %d\n",nvoxels);
  fflush(out_f);

  delete [] pt;
  delete [] pnt;
  delete [] vec;
  
  //then, begin cutting in next routine!
  //after cutting, main routine creates distinct nodes, boundaries to tetmesh, and then mesh!
  return;
  }
  
int split_tree::voxel_grad(List *nbrs0, List *nbrs1, int dir0, int dir1, double del[3], int i, int &changes3, double smn, double smx, Point pt1, Point pt2)
  {
  int flag = 0, k = 0, n0, n1, n2, n3;
  Point p0, p1, p2, p3;
  double deln[3];
  
  //both have 3 or 4 faces to my one, pick a direction to refine..will pick up 2-1 with next swipe
  if (nbrs0->max > 2 && nbrs1->max > 2)
    {
    //just an extra check
    if (voxels[i].children.max == 0)
      {
      refine_voxel(nvoxels, i, dir0+1);
      changes3++;
      flag = 1;
      }
    }
  else if (nbrs0->max > 2 && nbrs1->max == 2)
    {
    //only need to look at one neighbor
    p0 = voxels[nbrs1->list[0]].cornerpts[0];
    p1 = voxels[nbrs1->list[0]].cornerpts[1];
    for (k = 0; k < 3; k++)
      deln[k] = fabs(p1[k] - p0[k]);
            
    if (del[dir0] > 1.5*deln[dir0] && voxels[i].children.max == 0)
      {
      refine_voxel(nvoxels, i, dir0+1);
      changes3++;
      flag = 1;
      }
    else if (del[dir1] > 1.5*deln[dir1] && voxels[i].children.max == 0)
      {
      refine_voxel(nvoxels, i, dir1+1);
      changes3++;
      flag = 1;
      }
    }
  else if (nbrs0->max == 2 && nbrs1->max > 2)
    {
    //only need to look at one neighbor
    p0 = voxels[nbrs0->list[0]].cornerpts[0];
    p1 = voxels[nbrs0->list[0]].cornerpts[1];
    for (k = 0; k < 3; k++)
      deln[k] = fabs(p1[k] - p0[k]);
            
    if (del[dir0] > 1.5*deln[dir0] && voxels[i].children.max == 0)
      {
      refine_voxel(nvoxels, i, dir0+1);
      changes3++;
      flag = 1;
      }
    else if (del[dir1] > 1.5*deln[dir1] && voxels[i].children.max == 0)
      {
      refine_voxel(nvoxels, i, dir1+1);
      changes3++;
      flag = 1;
      }
    }
  else if (nbrs0->max == 2 && nbrs1->max == 2)
    {
    //we need both to agree, or else we create crossbar
    //due to the fact that we can have voxels bigger than the scope of our voxel, we cannot simply use deltas, or even cornerpts
    //so, we take the appropriate two pnts perturb two directions, and determine which yields 2 neighbors, and which only 1
    
    //perturb one side's point in the direction concerned with...use temp pts
    p0 = pt1;
    p1 = pt1;
    
    p0[dir0] += 0.25*smn;
    p1[dir0] -= 0.25*smn;
    
    //do nabor search
    n0 = nabor_search(p0,0,smn,smx);
    n1 = nabor_search(p1,0,smn,smx);
    
    //now, do other side
    p0 = pt2;
    p1 = pt2;
    
    p0[dir0] += 0.25*smn;
    p1[dir0] -= 0.25*smn;
    
    //do nabor search
    n2 = nabor_search(p0,0,smn,smx);
    n3 = nabor_search(p1,0,smn,smx);
    
    if (n0 != n1 && n2 != n3 && voxels[i].children.max == 0)
      {
      refine_voxel(nvoxels, i, dir0+1);
      changes3++;
      flag = 1;
      return(flag);
      }
    
    //perturb one side's point in the direction concerned with...use temp pts
    p0 = pt1;
    p1 = pt1;
    
    p0[dir1] += 0.25*smn;
    p1[dir1] -= 0.25*smn;
    
    //do nabor search
    n0 = nabor_search(p0,0,smn,smx);
    n1 = nabor_search(p1,0,smn,smx);
    
    //now, do other side
    p0 = pt2;
    p1 = pt2;
    
    p0[dir1] += 0.25*smn;
    p1[dir1] -= 0.25*smn;
    
    //do nabor search
    n2 = nabor_search(p0,0,smn,smx);
    n3 = nabor_search(p1,0,smn,smx);
    
    if (n0 != n1 && n2 != n3 && voxels[i].children.max == 0)
      {
      refine_voxel(nvoxels, i, dir1+1);
      changes3++;
      flag = 1;
      return(flag);
      }
      
    //if neither called, opp direction refined neighbors, and don't touch...flag still 0
    }
          
  return(flag);
  }

double sum_layer(int nl, double geomprog)
  {
  double length = 0.0;

  for (int i = 0; i < nl; i++)
    length += pow(geomprog,i);

  return(length);
  }
  
void split_tree::voxel_delete(geometry *geom, int in, double ctol, int f, int *vbd, int bd)
  {
  int n, m, flag = 0, flagv = 0;
  Point p0, p1, pn0, pn1, pn2, pn3, pg, pg0, pg1, pg2, pv;
  double begx, begy, begz, endx, endy, endz;
  //big enough to get edges of triangles and quads
  Point *pnt = new Point[4];
  Point *pt = new Point[4];
  Vector norm, v0, v1;
  double len = 0.0;
   
  //test for containment of pt/facet, intersection of edge with face in both directions which covers proximity/coplanar (uses fact that we have expanded voxel)
  //get extents
  p0 = voxels[in].cornerpts[0];
  p1 = voxels[in].cornerpts[1];
        
  //add in ctol extents
  p0 -= ctol;
  p1 += ctol;

  if (viscous)
    {
    flagv = 0;
    //check if bd is also viscous
    for (m = 0; m < nvbd && !flagv; m++)
      if (vbd[m] == bd)
        flagv = 1;
    if (flagv)
      {
      //calc normal from facet
      v0 = Vector(geom->g_vert[geom->g_facet[f].nodes[0]],geom->g_vert[geom->g_facet[f].nodes[1]]);
      v1 = Vector(geom->g_vert[geom->g_facet[f].nodes[0]],geom->g_vert[geom->g_facet[f].nodes[2]]);
      norm = v0 % v1;
      norm.normalize();
      //set viscous spacing
      len = sum_layer(viscous, factor);
      norm *= (vspace*len);
      }
    }
        
  for (n = 0; n < 3; n++)
    {
    if (viscous && flagv)
      {
      pg = geom->g_vert[geom->g_facet[f].nodes[n]];
      pg += Point(norm[0],norm[1],norm[2]);
      }
    else
      pg = geom->g_vert[geom->g_facet[f].nodes[n]];
          
    if (pg[0] < pn1[0] && pg[1] < pn1[1] && pg[2] < pn1[2] && pg[0] > pn0[0] && pg[1] > pn0[1] && pg[2] > pn0[2])
      {
      voxels[in].cut = 2;
      delete [] pt;
      delete [] pnt;
      return;
      }
    } 
              
  //set up extents to make new pts for next check, noting that all already expanded
  begx = p0[0];
  begy = p0[1];
  begz = p0[2];
  endx = p1[0];
  endy = p1[1];
  endz = p1[2];
          
  //now, check face intersections with facet edges (be sure to expand faces) and face edges with facets
  for (n = 0; n < 6; n++)
    {
    switch (n) 
      {
      case 0:
        pn0 = Point(begx,begy,begz);
        pn1 = Point(begx,endy,begz);
        pn2 = Point(endx,endy,begz);
        pn3 = Point(endx,begy,begz);
      break;
      case 1:
        pn0 = Point(begx,begy,begz);
        pn1 = Point(endx,begy,begz);
        pn2 = Point(endx,begy,endz);
        pn3 = Point(begx,begy,endz);
      break;
      case 2:
        pn0 = Point(endx,begy,begz);
        pn1 = Point(endx,endy,begz); 
        pn2 = Point(endx,endy,endz);
        pn3 = Point(endx,begy,endz);
      break;
      case 3:
        pn0 = Point(endx,endy,begz);
        pn1 = Point(begx,endy,begz);
        pn2 = Point(begx,endy,endz);
        pn3 = Point(endx,endy,endz);
      break;
      case 4:
        pn0 = Point(begx,endy,begz);
        pn1 = Point(begx,begy,begz);
        pn2 = Point(begx,begy,endz);
        pn3 = Point(begx,endy,endz);
      break;
      case 5:
        pn0 = Point(endx,endy,endz);
        pn1 = Point(begx,endy,endz);
        pn2 = Point(begx,begy,endz);
        pn3 = Point(endx,begy,endz);
      break;
      default:
      break;
      }
          
    //geom triangle
    if (viscous && flagv)
      {
      pg0 = geom->g_vert[geom->g_facet[f].nodes[0]];
      pg1 = geom->g_vert[geom->g_facet[f].nodes[1]];
      pg2 = geom->g_vert[geom->g_facet[f].nodes[2]];
      pg0 += Point(norm[0],norm[1],norm[2]);
      pg1 += Point(norm[0],norm[1],norm[2]);
      pg2 += Point(norm[0],norm[1],norm[2]);
      }
    else
      {
      pg0 = geom->g_vert[geom->g_facet[f].nodes[0]];
      pg1 = geom->g_vert[geom->g_facet[f].nodes[1]];
      pg2 = geom->g_vert[geom->g_facet[f].nodes[2]];
      } 
          
    //first set up triangle edges as line segments for loop
    pt[0] = pg0;
    pnt[0] = pg1;
    pt[1] = pg1;
    pnt[1] = pg2;
    pt[2] = pg2;
    pnt[2] = pg0;  
     
    //reset flag ..no need for external, since if it gets a flag here, it will return
    flag = 0;               
    for (m = 0; m < 3 && !flag; m++)
      {
      flag = line_facet_intersect(pt[m], pnt[m], pn0, pn1, pn2, pg, ctol);
      if (!flag)
        flag = line_facet_intersect(pt[m], pnt[m], pn2, pn3, pn0, pg, ctol);  
      }
            
    if (flag)
      {
      voxels[in].cut = 2;
      delete [] pt;
      delete [] pnt;
      return;
      }
      
    //first set up triangle edges as line segments for loop
    pt[0] = pn0;
    pnt[0] = pn1;
    pt[1] = pn1;
    pnt[1] = pn2;
    pt[2] = pn2;
    pnt[2] = pn3;
    pt[3] = pn3;
    pnt[3] = pn0;  
     
    //no need to reset, either 0 or returned (flag = 0)               
    for (m = 0; m < 4 && !flag; m++)
      {
      flag = line_facet_intersect(pt[m], pnt[m], pg0, pg1, pg2, pg, ctol);  
      }
            
    if (flag)
      {
      voxels[in].cut = 2;
      delete [] pt;
      delete [] pnt;
      return;
      }
    }

  delete [] pt;
  delete [] pnt;        
            
  return;        
  }
  
void split_tree::voxel_search(geometry *geom, double hi[3], double lo[3], int in, double ctol, int f, int *vbd, int bd)
  {
  int i;
  Point p0, p1;
  
  //first, set crnrpts
  p0 = voxels[in].cornerpts[0];
  p1 = voxels[in].cornerpts[1];
  
  //look to see if we're completely outside this voxel with ctol
  if (hi[0] < p0[0]-ctol || hi[1] < p0[1]-ctol || hi[2] < p0[2]-ctol || lo[0] > p1[0]+ctol || lo[1] > p1[1]+ctol || lo[2] > p1[2]+ctol)
    {
    return;
    }
  else
    {
    if (voxels[in].children.max > 0)
      {
      for (i = 0; i < voxels[in].children.max; i++)
        {
        voxel_search(geom,hi,lo,voxels[in].children.list[i],ctol,f, vbd, bd);
        }
      }
    else
      {
      if (f == -1)
        voxels[in].cut = 2;
      else
        voxel_delete(geom,in,ctol,f,vbd, bd);
      }
    }

  return;
  }

void split_tree::voxel_search_tree(Point tpt, int in, double ctol, List *out)
  {
  int i;
  Point p0, p1;
  
  //first, set crnrpts
  p0 = voxels[in].cornerpts[0];
  p1 = voxels[in].cornerpts[1];
  
  //look to see if we're completely outside this voxel with ctol
  if (tpt[0] < p0[0]-ctol || tpt[1] < p0[1]-ctol || tpt[2] < p0[2]-ctol || tpt[0] > p1[0]+ctol || tpt[1] > p1[1]+ctol || tpt[2] > p1[2]+ctol)
    {
    return;
    }
  else
    {
    if (voxels[in].children.max > 0)
      {
      for (i = 0; i < voxels[in].children.max; i++)
        {
        voxel_search_tree(tpt,voxels[in].children.list[i],ctol,out);
        }
      }
    else
      {
      if (voxels[in].cut == 0)
        out->Check_List(in);
      }
    }

  return;
  }
  
void split_tree::cut_voxels(geometry *geom, double ctol, double smn, double smx, int *boxcut, int *vbd)
  {
  int i, j, k, n, m, l, n0, n1, n2, flag = 0, flagv = 0;
  Point p0, p1, cent, pv;
  double hi[3],lo[3];
  Vector norm, v0, v1;
  time_t tm;
  char *t_char;
  double len = 0.0;

  time(&tm);
  t_char = ctime(&tm);
  fprintf(out_f,"\nBeginning voxel cutting at %s.\n",t_char);
  fflush(out_f);
  
  //if cutting out a full box of geom
  if (evencut)
    {
    for (i = 0; i < bb; i++)
      {
      //reset hi, lo
      for (m = 0; m < 3; m++)
        {
        hi[m] = -1.0e+20;
        lo[m] = 1.0e+20;
        }
      for (j=geom->n_begin[boxcut[i]]; j < geom->n_end[boxcut[i]]; j++)
        {
        //first, determine geom extents
        n0 = geom->g_facet[j].nodes[0];
        n1 = geom->g_facet[j].nodes[1];
        n2 = geom->g_facet[j].nodes[2];

        if (viscous)
          {
          flagv = 0;
          //check if boxcut bd is also viscous
          for (m = 0; m < nvbd && !flagv; m++)
            if (vbd[m] == boxcut[i])
              flagv = 1;
          if (flagv)
            {
            //calc normal from facet
            v0 = Vector(geom->g_vert[n0],geom->g_vert[n1]);
            v1 = Vector(geom->g_vert[n0],geom->g_vert[n2]);
            norm = v0 % v1;
            norm.normalize();
            //set viscous spacing
            len = sum_layer(viscous, factor);
            norm *= (vspace*len);
            }
          }

        //now, set up facet extents for use in determining where to go in children
        for (n = 0; n < 3; n++)
          {
          switch (n) 
            {
            case 0:
              k = n0;
            break;
            case 1:
              k = n1;
            break;
            case 2:
              k = n2;
            break;
            default:
              fprintf(stderr,"\nFacet %d has more than three sides.  Exiting....\n",j);
              fflush(stderr);
              //MPI_Abort(MPI_COMM_WORLD,0);
              exit(0);
            break;
            }
          for (m = 0; m < 3; m++)
            {
            if (viscous && flagv)
              {
              pv = Point(geom->g_vert[k][0],geom->g_vert[k][1],geom->g_vert[k][2]);
              pv += Point(norm[0],norm[1],norm[2]);

              hi[m] = MAX(hi[m],pv[m]);
              lo[m] = MIN(lo[m],pv[m]);
              }
            else
              {
              hi[m] = MAX(hi[m],geom->g_vert[k][m]);
              lo[m] = MIN(lo[m],geom->g_vert[k][m]);
              }
            }
          }
        }
      //now, recursively find voxels to cut (to avoid making a huge list)
      voxel_search(geom, hi, lo, 0, ctol, -1, vbd, boxcut[i]);
      }
    }

  //need to loop through looking if in extents, then try intersection both ways and will dump voxels too close to geom (takes care of too close to each other, as well)
  for (i=1; i <= geom->ngb; i++)
    {
    //check for boxcut
    flag = 0;
    for (j = 0; j < bb && !flag; j++)
      if (boxcut[j] == i)
        flag = 1;

    if (flag)
      continue;

    for (j=geom->n_begin[i]; j < geom->n_end[i]; j++)
      {
      //first, determine geom extents
      n0 = geom->g_facet[j].nodes[0];
      n1 = geom->g_facet[j].nodes[1];
      n2 = geom->g_facet[j].nodes[2];

      if (viscous)
        {
        flagv = 0;
        //check if bd is also viscous
        for (m = 0; m < nvbd && !flagv; m++)
          if (vbd[m] == i)
            flagv = 1;
        if (flagv)
          {
          //calc normal from facet
          v0 = Vector(geom->g_vert[n0],geom->g_vert[n1]);
          v1 = Vector(geom->g_vert[n0],geom->g_vert[n2]);
          norm = v0 % v1;
          norm.normalize();
          //set viscous spacing
          len = sum_layer(viscous, factor);
          norm *= (vspace*len);
          }
        }
      
      //reset hi, lo
      for (m = 0; m < 3; m++)
        {
        hi[m] = -1.0e+20;
        lo[m] = 1.0e+20;
        }
          
      //now, set up facet extents for use in determining where to go in children
      for (n = 0; n < 3; n++)
        {
        switch (n) 
          {
          case 0:
            k = n0;
          break;
          case 1:
            k = n1;
          break;
          case 2:
            k = n2;
          break;
          default:
            fprintf(stderr,"\nFacet %d has more than three sides.  Exiting....\n",j);
            fflush(stderr);
            //MPI_Abort(MPI_COMM_WORLD,0);
            exit(0);
          break;
          }
        for (m = 0; m < 3; m++)
          {
          if (viscous && flagv)
            {
            pv = Point(geom->g_vert[k][0],geom->g_vert[k][1],geom->g_vert[k][2]);
            pv += Point(norm[0],norm[1],norm[2]);

            hi[m] = MAX(hi[m],pv[m]);
            lo[m] = MIN(lo[m],pv[m]);
            }
          else
            {
            hi[m] = MAX(hi[m],geom->g_vert[k][m]);
            lo[m] = MIN(lo[m],geom->g_vert[k][m]);
            }
          }
        }
        
      //now, recursively find voxels to cut (to avoid making a huge list)
      //once found in search, will call voxel_delete
      voxel_search(geom, hi, lo, 0, ctol, j, vbd, i);
      }
    }

  time(&tm);
  t_char = ctime(&tm);
  fprintf(out_f,"\nFinished voxel cutting at %s.\n",t_char);
  fflush(out_f);

  flag = 0;
  for (i = 0; i < nvoxels; i++)
    if (voxels[i].children.max == 0 && voxels[i].cut == 2)
      flag++;
  
  fprintf(out_f,"\nFinished cutting with %d voxels marked as cut.  Beginning flood fill.\n",flag);
  fflush(out_f);
  
  //determine how many voxels need to be marked for flood fill to be complete
  int flag1 = 0, flag2 = 0;
  for (i = 0; i < nvoxels; i++)
    {
    if (voxels[i].children.max == 0)
      {
      flag2++;
      if (voxels[i].cut != 2)
        flag1++;
      }
    }
  
  fprintf(out_f,"\nNeed to flood fill %d out of %d voxels at finest level.\n",flag1,flag2);
  fflush(out_f);
  
  time(&tm);
  t_char = ctime(&tm);
  fprintf(out_f,"\nBeginning flood-fill at %s.\n",t_char);
  fflush(out_f);

  //finally, flood fill and use a do loop to be sure all gotten
  //cannot use a flood fill recursion since walks through whole mesh if all directions taken and stack gets too big
  //on even bigger meshes, will run out even if all directions not taken
  //use an "infinite" loop instead to save stack overhead
  int changes = 0, intchanges = 0;
  //int tempn = 0;
  int pass = 0;
  List noreturn; //will hold voxels whose in/out is indeterminate so filled by neighbors
  noreturn.construct();
  Point *pt = new Point[6]; 
  Point *pnt = new Point[6];
  Vector *vec = new Vector[6];
  int nbr[6];
  double begx, begy, begz, endx, endy, endz, halfx, halfy, halfz;
  //not doing corners so no danger of jumping across partition
  double tol = 0.25*smn;
  
  if (!altflood)
  {
  do
    {
    pass++;
    flag = flag2 = 0; 
    for (i = 0; i < nvoxels && !flag; i++)
      {
      if (voxels[i].children.max == 0 && voxels[i].cut == 0 && !noreturn.Is_In_List(i))
        {
        n = i;
        flag = 1;
        }
      }
        
    //now, get this voxel's in out status and get as many neighbors as you can
    //get crnerpts
    p0 = voxels[n].cornerpts[0];
    p1 = voxels[n].cornerpts[1];
  
    //get centroid
    cent = p0 + p1;
    cent /= 2.0;
      
    //now, look if voxel is in or out, find neighbors, flood fill
    //will give in (1), out (-1), or unknown (left 0)
    if (voxels[n].cut == 0)
      {
      flag2++;
      voxels[n].cut = geom->In_Out(cent, 1.0e-15);
      }
      
    if (voxels[n].cut == 0)
      {
      noreturn.Add_To_List(n); //we got no information from this voxel
      //fprintf(out_f,"\nVoxel %d gave no info.\n",n);
      //fflush(out_f);
      }
    else if (flag2)
      {
      changes++;
      //reset internal changes to go through loop
      intchanges = 1;
      do
        {
        //reset internal changes in loop
        intchanges = 0;
        
        //set current voxel to get neighbors
        i = n;
        
        //get crnerpts
        p0 = voxels[i].cornerpts[0];
        p1 = voxels[i].cornerpts[1];
  
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
        
        flag = 0; //reset to jump out when we find a neighbor
        //for (j = 0; j < 6; j++)
        for (j = 0; j < 6 && !flag; j++)
          {
          vec[j] = Vector(cent,pt[j]); //get direction
          vec[j].normalize(); //normalize
          vec[j] *= tol; //get tol as length
          pnt[j] = pt[j] + Point(vec[j][0],vec[j][1],vec[j][2]);
          nbr[j] = nabor_search(pnt[j],0,smn,smx);
      
          if (nbr[j] < 0 || voxels[nbr[j]].cut != 0 || voxels[nbr[j]].children.max > 0)
            continue; //no need to mess with non-existant or already marked as cut, in, out
        
          //here, the assumption is that there will be a full voxel band between two domains
          if (voxels[nbr[j]].cut == 0)
            {  
            voxels[nbr[j]].cut = voxels[i].cut; //set as voxel in question
            changes++;
            intchanges++;
            flag++;
            n = nbr[j];
            //tempn = nbr[j]; //holds most recent legit neighbor
            }
          else if (voxels[nbr[j]].cut != 0 && voxels[nbr[j]].cut != voxels[i].cut)
            {
            fprintf(out_f,"\nVoxel %d disagrees with voxel %d about in/out status.  Exiting....\n",nbr[j],i);
            voxels[nbr[j]].print(out_f);
            voxels[i].print(out_f);
            fflush(out_f);
            exit(0);
            }
          }
          
        //if (tempn >= 0)
          //n = tempn;  
          
        }while (intchanges > 0);
      }

    fprintf(out_f,"\nFlood fill pass = %d, changes = %d.\n",pass,changes);
    fflush(out_f);
    }while (changes < flag1);
  }
  else
  {
  time(&tm);
  t_char = ctime(&tm);
  fprintf(out_f,"\nStarting nabor table at %s.\n",t_char);
  fflush(out_f);
  //make a neighbor table and flood fill that way
  int **nbrtbl = 0;
  nbrtbl = new int*[nvoxels];
  for (i = 0; i < nvoxels; i++)
    {
    nbrtbl[i] = new int [6];
    for (j = 0; j < 6; j++)
      nbrtbl[i][j] = -1;
    }

  for (i = 0; i < nvoxels; i++)
    {
    if (voxels[i].children.max > 0 || voxels[i].cut == 2)
      continue;
 
    //get crnerpts
    p0 = voxels[i].cornerpts[0];
    p1 = voxels[i].cornerpts[1];
  
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
        
    for (j = 0; j < 6; j++)
      {
      vec[j] = Vector(cent,pt[j]); //get direction
      vec[j].normalize(); //normalize
      vec[j] *= tol; //get tol as length
      pnt[j] = pt[j] + Point(vec[j][0],vec[j][1],vec[j][2]);
      nbrtbl[i][j] = nabor_search(pnt[j],0,smn,smx);
      }
    }

  time(&tm);
  t_char = ctime(&tm);
  fprintf(out_f,"\nFinished nabor table at %s.\n",t_char);
  fflush(out_f);

  if (floodnum)
    {
    //create a grid, pass pts down tree, ray cast
    Point *tpts = new Point[floodnum*floodnum*floodnum];
    for (i = 0; i < floodnum*floodnum*floodnum; i++)
      tpts[i] = Point(0.0,0.0,0.0);

    //find extents
    p0 = geom->min_point();
    p1 = geom->max_point();

    //find deltas
    double deltax = fabs(p1[0] - p0[0])/(floodnum-1);
    double deltay = fabs(p1[1] - p0[1])/(floodnum-1);
    double deltaz = fabs(p1[2] - p0[2])/(floodnum-1);

    //now, create tpts
    l = 0;
    for (i = 0; i < floodnum; i++)
      {
      for (j = 0; j < floodnum; j++)
        {
        for (k = 0; k < floodnum; k++)
          {
          tpts[l] = Point(p0[0]+i*deltax, p0[1]+j*deltay, p0[2]+k*deltaz);
          l++;
          }
        }
      }

    List *out = new List();
    time(&tm);
    t_char = ctime(&tm);
    fprintf(out_f,"\nStarting tree traversal at %s.\n",t_char);
    fflush(out_f);
    //now, pass each pt down tree, find and check voxel, do in out test
    for (i = 0; i < floodnum*floodnum*floodnum; i++)
      voxel_search_tree(tpts[i], 0, ctol, out);

    //debug
    //out->print(out_f);

    time(&tm);
    t_char = ctime(&tm);
    fprintf(out_f,"\nFinished tree traversal at %s.\n",t_char);
    fflush(out_f);

    for (i = 0; i < out->max; i++)
      {
      p0 = voxels[out->list[i]].cornerpts[0];
      p1 = voxels[out->list[i]].cornerpts[1];
      cent = p0 + p1;
      cent /= 2.0;
      int t = geom->In_Out(cent, 1.0e-15);
      if (t != 0 && voxels[out->list[i]].cut == 0)
        {
        voxels[out->list[i]].cut = t;
        //fprintf(out_f,"\nMarked voxel %d\n",out->list[i]);
        //fflush(out_f);
        changes++;
        }
      }

    time(&tm);
    t_char = ctime(&tm);
    fprintf(out_f,"\nFinished %d ray traces at %s, yielding %d changes.\n",out->max,t_char,changes);
    fflush(out_f);

    delete out;
    delete [] tpts;
    }

  //now, get voxel's in/out status, and use nabors to flood fill
  do
    {
    pass++;
    time(&tm);
    t_char = ctime(&tm);
    fprintf(out_f,"\nPass %d started at %s.\n",pass,t_char);
    fflush(out_f);
    flag = flag2 = 0; 
    for (i = 0; i < nvoxels && !flag; i++)
      {
      if (voxels[i].children.max == 0 && voxels[i].cut == 0 && !noreturn.Is_In_List(i))
        {
        n = i;
        flag = 1;
        }
      }

    if (!flag)
      {
      fprintf(out_f,"\nYou've run out of flood-fillable voxels.  Exiting....\n");
      fprintf(out_f,"\nnoreturn:\n");
      noreturn.print(out_f);
      fflush(out_f);
      exit(0);
      }
        
    //now, get this voxel's in out status and get as many neighbors as you can
    //get crnerpts
    p0 = voxels[n].cornerpts[0];
    p1 = voxels[n].cornerpts[1];
  
    //get centroid
    cent = p0 + p1;
    cent /= 2.0;
      
    //now, look if voxel is in or out, find neighbors, flood fill
    //will give in (1), out (-1), or unknown (left 0)
    if (voxels[n].cut == 0)
      {
      flag2++;
      voxels[n].cut = geom->In_Out(cent, 1.0e-15);
      }
      
    if (voxels[n].cut == 0)
      {
      noreturn.Add_To_List(n); //we got no information from this voxel
      fprintf(out_f,"\nVoxel %d gave no info...no return max = %d.\n",n,noreturn.max);
      fflush(out_f);
      }
    else if (flag2)
      {
      changes++; //we demarcated one voxel
      time(&tm);
      t_char = ctime(&tm);
      fprintf(out_f,"\nInner loops for pass %d started at %s.\n",pass,t_char);
      fflush(out_f);
      //now, inform nabors
      do
        {
        intchanges = 0; //reset so when cannot fill any more, go find other ray trace voxel
        for (i = 0; i < nvoxels; i++)
          {
          if (voxels[i].children.max > 0 || voxels[i].cut == 0 || voxels[i].cut == 2)
            continue; //not at finest level, undetermined, or cut...move on
          //otherwise, we have a marked voxel and need to find its nabors
          for (j = 0; j < 6; j++)
            {
            if (nbrtbl[i][j] < 0 || voxels[nbrtbl[i][j]].children.max > 0 || voxels[nbrtbl[i][j]].cut == voxels[i].cut || voxels[nbrtbl[i][j]].cut == 2)
              continue; //no nabor, has kids, already set, or cut
            else if (voxels[nbrtbl[i][j]].cut != 0 && voxels[nbrtbl[i][j]].cut != voxels[i].cut)
              {
              //check for disagreement
              fprintf(out_f,"\nVoxel %d is marked %d and its nabor %d is marked %d.  Exiting....\n",i,voxels[i].cut,nbrtbl[i][j],voxels[nbrtbl[i][j]].cut);
              fflush(out_f);
              exit(0);
              }
            else
              {
              voxels[nbrtbl[i][j]].cut = voxels[i].cut;
              intchanges++;
              changes++;
              }
            }
          }
        } while (intchanges > 0);

      time(&tm);
      t_char = ctime(&tm);
      fprintf(out_f,"\nInner loops for pass %d ended at %s.\n",pass,t_char);
      fflush(out_f);
      }
    time(&tm);
    t_char = ctime(&tm);
    fprintf(out_f,"\nFlood fill pass = %d, changes = %d, ended at %s.\n",pass,changes,t_char);
    fflush(out_f);
    }while (changes < flag1);

  for (i = 0; i < nvoxels; i++)
    delete [] nbrtbl[i];
  delete [] nbrtbl;
  }

  time(&tm);
  t_char = ctime(&tm);
  fprintf(out_f,"\nFinished flood-fill at %s.\n",t_char);
  fflush(out_f);
    
  //DEBUG!  
  flag1 = flag2 = 0;
  for (i = 0; i < nvoxels; i++)
    if (voxels[i].cut == 0 && voxels[i].children.max == 0)
      {
      flag1++;
      if (noreturn.Is_In_List(i))
        flag2++;
      }
      
  noreturn.destruct();

  fprintf(out_f,"\n%d voxels at finest level unmarked, with %d in the noreturn list (max = %d).\n",flag1,flag2,noreturn.max);
  fflush(out_f);
  
  flag1 = 0;
  for (i = 0; i < nvoxels; i++)
    if (voxels[i].cut != 0 && voxels[i].children.max > 0)
      flag1++;

  fprintf(out_f,"\n%d voxels NOT at finest level marked.\n",flag1);
  fflush(out_f);
    
  delete [] vec;
  delete [] pt;
  delete [] pnt;
  
  fprintf(out_f,"\nFinshed cutting and flood fill.\n");
  fflush(out_f);
  
  return;
  }
  
  int split_tree::eval_qualityCB(int v, int nbr, double del[3], int type[2], int &changes)
    {
    int flag = 0;
    double deln[2];

    //find neighbor deltas
    deln[type[0]-1] = fabs(voxels[nbr].cornerpts[1][type[0]-1] - voxels[nbr].cornerpts[0][type[0]-1]);
    deln[type[1]-1] = fabs(voxels[nbr].cornerpts[1][type[1]-1] - voxels[nbr].cornerpts[0][type[1]-1]);
    
    //now, perform tests with proper neighbors...look for 4-1, then crossbar
    if (deln[type[0]-1] > 1.5*del[type[0]-1] && 1.5*deln[type[1]-1] < del[type[1]-1] && voxels[nbr].children.max == 0)
      {
      refine_voxel(nvoxels, nbr, type[0]);
      if (voxels[v].children.max == 0)
        {
        refine_voxel(nvoxels, v, type[1]);
        }
      flag = 1;
      changes++;
      }
    else if (deln[type[1]-1] > 1.5*del[type[1]-1] && 1.5*deln[type[0]-1] < del[type[0]-1] && voxels[nbr].children.max == 0)
      {
      refine_voxel(nvoxels, nbr, type[1]);
      if (voxels[v].children.max == 0)
        {
        refine_voxel(nvoxels, v, type[0]);
        }
      flag = 1;
      changes++;
      }
    return(flag);
    }
    
  void split_tree::eval_quality41(int nbr, double del[3], int type[2])
    {
    int i;
    double deln[3];

    //find neighbor deltas
    for (i = 0; i < 3; i++)
      deln[i] = fabs(voxels[nbr].cornerpts[1][i] - voxels[nbr].cornerpts[0][i]);
    
    if (deln[type[0]-1] > 3.0*del[type[0]-1] && voxels[nbr].children.max == 0) 
      {
      if (voxels[nbr].cut == 0 || voxels[nbr].cut > type[0])
        voxels[nbr].cut = type[0];
      }
    else if (deln[type[1]-1] > 3.0*del[type[1]-1] && voxels[nbr].children.max == 0)
      { 
      if (voxels[nbr].cut == 0 || voxels[nbr].cut > type[1])
        voxels[nbr].cut = type[1];
      }
    return;
    }
