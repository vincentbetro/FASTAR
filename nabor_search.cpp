#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "split_tree.h"

#define ABS(x) ((x) >= 0 ? (x) : -(x))
#define SIGN(a,b) ((b) < 0.0 ? -ABS(a) : ABS(a))

//global variables
extern int my_rank; /*rank of process*/
extern int num_procs; /*number of processes*/
extern int adapt; //using adaptation or not
extern FILE *in_f, *jou_f, *out_f, *debug_f; /*input output journal files global*/

//returns nabor or -1 if none
int split_tree::nabor_search(Point pt, int in, double smn, double smx)
  {
  int s, out, i, j, m;
  Point p0, p1, cent;
  
  //this needs to be min over max times min, since this is the least amount of smn 
  //that could be contributed to a cardinal direction in a non cardinal unit vector of length smn
  //should mult by 0.25 since smn can be 0.5 of smn, but figured smaller was better!
  double tol =(smn/smx)*smn*0.01;

  //start at root
  p0 = voxels[in].cornerpts[0];
  p1 = voxels[in].cornerpts[1];

  //check to be sure not outside of root, if so, return -1
  //okay to use cornerpts since off root, which has actual physical coords stored
  if (pt[0] > voxels[0].cornerpts[1][0] || pt[1] > voxels[0].cornerpts[1][1] || pt[2] > voxels[0].cornerpts[1][2] || pt[0] < voxels[0].cornerpts[0][0] || pt[1] < voxels[0].cornerpts[0][1] || pt[2] < voxels[0].cornerpts[0][2])
    return(-1);
        
  //compute centroid
  cent = p0 + p1;
  cent /= 2.0;
  
  //now, determine next direction to go based on split
  s = voxels[in].split;
  
  switch (s) 
    {
    case 0:
      //we've found the neighbor
      return(in);
    break;
    case 1:
      i = voxels[in].children.list[0];
      j = voxels[in].children.list[1];
      if (pt[0] < (cent[0]-tol))
        out = i;
      else if (pt[0] > (cent[0]+tol))
        out = j;
      else
        return (in);
    break;
    case 2:
      i = voxels[in].children.list[0];
      j = voxels[in].children.list[1];
      if (pt[1] < (cent[1]-tol))
        out = i;
      else if (pt[1] > (cent[1]+tol))
        out = j;
      else
        return (in);
    break;
    case 3:
      i = voxels[in].children.list[0];
      j = voxels[in].children.list[1];
      if (pt[2] < (cent[2]-tol))
        out = i;
      else if (pt[2] > (cent[2]+tol))
        out = j;
      else
        return (in);
    break;
    default:
      fprintf(stderr,"\nYou have an undefined split for voxel %d.  Exiting....\n",in);
      fflush(stderr);
      exit(0);
    break;
    }
  
  m = nabor_search(pt,out,smn,smx);

  return(m);
  }
  
void split_tree::set_nabors_nodes(int k, int nbr, double tol, Point pt)
  {
  int n, flag = 0;
  Point np0, np1, ncent;
  Point *npt = new Point[26];
  
  //get crnerpts
  np0 = voxels[nbr].cornerpts[0];
  np1 = voxels[nbr].cornerpts[1];
  
  //get centroid
  ncent = np0 + np1;
  ncent /= 2.0;
  
  nabor_nodes(npt, np0, np1);
      
  //now, look at approprate neighbors and determine if they have the node in question
  //be sure if neighbor doesn't have a slot, we kick out!!!!
  for (n = 0; n < 26 && !flag; n++)
    {
    if (npt[n][0] < pt[0]+tol && npt[n][0] > pt[0]-tol && npt[n][1] < pt[1]+tol && npt[n][1] > pt[1]-tol && npt[n][2] < pt[2]+tol && npt[n][2] > pt[2]-tol)
      {
      flag = 1;
      if (voxels[nbr].nodes[n] < 0)
        {
        voxels[nbr].nodes[n] = k; 
        }
      else if (voxels[nbr].nodes[n] >= 0 && voxels[nbr].nodes[n] != k)
        {
        fprintf(stderr,"\n2)You have found a neighbor (%d) (%lf, %lf, %lf) (%lf, %lf, %lf) of a voxel that already has the given neighboring point (%lf, %lf, %lf) at index %d.  Exiting....\n",nbr,voxels[nbr].cornerpts[0][0],voxels[nbr].cornerpts[0][1],voxels[nbr].cornerpts[0][2],voxels[nbr].cornerpts[1][0],voxels[nbr].cornerpts[1][1],voxels[nbr].cornerpts[1][2],pt[0],pt[1],pt[2],n);
        voxels[nbr].print(stderr);
        fflush(stderr);
        exit(0);
        }
      }
    }

  delete [] npt;

  return;
  }
  
int split_tree::find_nabors_nodes(int nbr, double tol, Point pt)
  {
  int n, flag = 0;
  Point np0, np1, ncent;
  Point *npt = new Point[26];
  
  //get crnerpts
  np0 = voxels[nbr].cornerpts[0];
  np1 = voxels[nbr].cornerpts[1];
  
  //get centroid
  ncent = np0 + np1;
  ncent /= 2.0;
  
  nabor_nodes(npt, np0, np1);
      
  //now, look at approprate neighbors and determine if they have the node in question
  //be sure if neighbor doesn't have a slot, we kick out!!!!
  for (n = 0; n < 26 && !flag; n++)
    {
    if (npt[n][0] < pt[0]+tol && npt[n][0] > pt[0]-tol && npt[n][1] < pt[1]+tol && npt[n][1] > pt[1]-tol && npt[n][2] < pt[2]+tol && npt[n][2] > pt[2]-tol)
      {
      flag = 1;
      if (voxels[nbr].nodes[n] >= 0)
        {
        delete [] npt;
        return(voxels[nbr].nodes[n]); 
        }
      }
    }
  if (!flag)
    { 
    fprintf(stderr,"\n1)You have found a neighbor (%d) (%lf, %lf, %lf) (%lf, %lf, %lf) of a voxel that cannot contain the given neighboring point (%lf, %lf, %lf).  Exiting....\n",nbr,voxels[nbr].cornerpts[0][0],voxels[nbr].cornerpts[0][1],voxels[nbr].cornerpts[0][2],voxels[nbr].cornerpts[1][0],voxels[nbr].cornerpts[1][1],voxels[nbr].cornerpts[1][2],pt[0],pt[1],pt[2]);
    voxels[nbr].print(stderr);
    fflush(stderr);
    exit(0);
    }

  delete [] npt;

  return (-1);
  }

void split_tree::set_nodes(int v, double smn, double smx, int &num_nodes, int *nums, Point *pts)
  {
  int j, k, n, flag = 0;
  Point p0, p1, cent;
  double tol;
  Point *pt = new Point[26]; 
  Point *pnt = new Point[26];
  Vector *vec = new Vector[26];
  int nbr[26];
  int nbrtest[7];
    
  //we need to check all directions
  //get crnerpts
  p0 = voxels[v].cornerpts[0];
  p1 = voxels[v].cornerpts[1];
  
  //get centroid
  cent = p0 + p1;
  cent /= 2.0;
  
  nabor_nodes(pt, p0, p1);
  
  //set up tolerance length
  tol = 0.25*smn;
      
  //find all neigbors    
  for (j = 0; j < 26; j++)
    {
    vec[j] = Vector(cent,pt[j]); //get direction
    vec[j].normalize(); //normalize
    vec[j] *= tol; //get tol as length
    pnt[j] = pt[j] + Point(vec[j][0],vec[j][1],vec[j][2]);
    nbr[j] = nabor_search(pnt[j],0,smn,smx);
    }
  
  //now, look through each cornerpt and all neighbors that may have
  for (j = 0; j < 8; j++)
    {
    if (voxels[v].nodes[j] >= 0)
      continue; //don't look if already found
    
    flag = 0;
    switch (j) 
      {
      case 0:
        nbrtest[0] = nbr[0];
        nbrtest[1] = nbr[8];
        nbrtest[2] = nbr[11];
        nbrtest[3] = nbr[12];
        nbrtest[4] = nbr[20];
        nbrtest[5] = nbr[21];
        nbrtest[6] = nbr[24];
      break;
      case 1:
        nbrtest[0] = nbr[1];
        nbrtest[1] = nbr[8];
        nbrtest[2] = nbr[9];
        nbrtest[3] = nbr[13];
        nbrtest[4] = nbr[20];
        nbrtest[5] = nbr[21];
        nbrtest[6] = nbr[22];
      break;
      case 2:
        nbrtest[0] = nbr[2];
        nbrtest[1] = nbr[9];
        nbrtest[2] = nbr[10];
        nbrtest[3] = nbr[14];
        nbrtest[4] = nbr[20];
        nbrtest[5] = nbr[22];
        nbrtest[6] = nbr[23];
      break;
      case 3:
        nbrtest[0] = nbr[3];
        nbrtest[1] = nbr[10];
        nbrtest[2] = nbr[11];
        nbrtest[3] = nbr[15];
        nbrtest[4] = nbr[20];
        nbrtest[5] = nbr[23];
        nbrtest[6] = nbr[24];
      break;
      case 4:
        nbrtest[0] = nbr[4];
        nbrtest[1] = nbr[12];
        nbrtest[2] = nbr[16];
        nbrtest[3] = nbr[19];
        nbrtest[4] = nbr[25];
        nbrtest[5] = nbr[21];
        nbrtest[6] = nbr[24];
      break;
      case 5:
        nbrtest[0] = nbr[5];
        nbrtest[1] = nbr[13];
        nbrtest[2] = nbr[16];
        nbrtest[3] = nbr[17];
        nbrtest[4] = nbr[25];
        nbrtest[5] = nbr[21];
        nbrtest[6] = nbr[22];
      break;
      case 6:
        nbrtest[0] = nbr[6];
        nbrtest[1] = nbr[14];
        nbrtest[2] = nbr[17];
        nbrtest[3] = nbr[18];
        nbrtest[4] = nbr[25];
        nbrtest[5] = nbr[22];
        nbrtest[6] = nbr[23];
      break;
      case 7:
        nbrtest[0] = nbr[7];
        nbrtest[1] = nbr[15];
        nbrtest[2] = nbr[18];
        nbrtest[3] = nbr[19];
        nbrtest[4] = nbr[25];
        nbrtest[5] = nbr[23];
        nbrtest[6] = nbr[24];
      break;
      default:
      break;
      }
    for (k = 0; k < 7 && !flag; k++)
      {
      //only look at neighbor if it is at finest level..else, let finest level determine all pts
      if (nbrtest[k] < 0 || voxels[nbrtest[k]].children.max > 0 || voxels[nbrtest[k]].cut <= 0)
        continue; //will get this bigger voxel sorted by littler ones
        
      n = find_nabors_nodes(nbrtest[k], tol, pt[j]);

      if (n >= 0)
        {
        flag = 1;
        }
      }
      
    //if found, take number and give to other neighbors
    if (flag)
      {
      voxels[v].nodes[j] = n;
      
      for (k = 0; k < 7; k++)
        {
        //only look at neighbor if it is at finest level..else, let finest level determine all pts
        if (nbrtest[k] < 0 || voxels[nbrtest[k]].children.max > 0 || voxels[nbrtest[k]].cut <= 0)
          continue; //will get this bigger voxel sorted by littler ones
        
        set_nabors_nodes(n, nbrtest[k], tol, pt[j]);      
        }
      }
    else
      {
      //if not, name node, let neighbors know
      //set voxel node number
      voxels[v].nodes[j] = num_nodes;
      //set up node number in nums
      nums[j] = num_nodes;
      //set up physical coord in pts
      pts[j] = pt[j];
      //increment num_nodes for return
      num_nodes++;
      
      //pass to neighbors
      for (k = 0; k < 7; k++)
        {
        //only look at neighbor if it is at finest level..else, let finest level determine all pts
        if (nbrtest[k] < 0 || voxels[nbrtest[k]].children.max > 0 || voxels[nbrtest[k]].cut <= 0)
          continue; //will get this bigger voxel sorted by littler ones
        
        set_nabors_nodes(nums[j], nbrtest[k], tol, pt[j]);      
        }
      }
    }

  delete [] pt;
  delete [] pnt;
  delete [] vec;

  return;                
  }
  
void split_tree::nabor_nodes(Point *pt, Point p0, Point p1)
  {
  double begx, begy, begz, endx, endy, endz, halfx, halfy, halfz;
  
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
  pt[0] = Point(begx,begy,begz);
  pt[1] = Point(endx,begy,begz);
  pt[2] = Point(endx,endy,begz);
  pt[3] = Point(begx,endy,begz);
  pt[4] = Point(begx,begy,endz);
  pt[5] = Point(endx,begy,endz);
  pt[6] = Point(endx,endy,endz);
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
  pt[20] = Point(halfx,halfy,begz);
  pt[21] = Point(halfx,begy,halfz);
  pt[22] = Point(endx,halfy,halfz);
  pt[23] = Point(halfx,endy,halfz);
  pt[24] = Point(begx,halfy,halfz); 
  pt[25] = Point(halfx,halfy,endz);
  
  return;
  }
  
