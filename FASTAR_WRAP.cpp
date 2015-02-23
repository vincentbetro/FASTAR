#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include "geometry.h"
#include "FASTAR.h"
#include "Vector.h"
#include "journal.h"
#include "Spacing_Field.h"
#define MIN(x,y) ((x) <= (y) ? (x) : (y))
#define MAX(x,y) ((x) >= (y) ? (x) : (y))

//global variables
int my_rank; /*rank of process*/
int num_procs; /*number of processes*/
int adapt; //using adaptation or not
int restart; //restarting or not
int ugeom; //spacing file only (1) or both spacing and geom used
int ftype; //0 for UCD, 1 for STL
int evencut; //if we want nice box around geometry
int bb; //num bd to box cut
int gentet; //use tetgen
int viscous; //insert layers
double factor; //geom prog factor
double vspace; //viscous spacing
int nvbd; //number viscous bd
int highquality; //allows desirable constraints to be enforced
int altflood; //type of flood fill
int floodnum; //cube root of num pts to do ray trace
int numregions; //gives number of regions to keep from tetgen
int overset; //overset or not
FILE *in_f, *jou_f, *out_f, *debug_f; /*input output journal files global*/

int main(int argcs, char* pArgs[])
{
  //declare variables (if not used as a library)
  //int source; //rank of sender
  //int dest; //rank of receiver
  //int tag = 0; //tag for messages
  //MPI_Status status; //return status for receive
  int i, nf;
  //int digits; //counters and input params
  const int bdim = 132; //buffer dim
  //char extension[bdim]; //file extension to allow padded digits
  char** fnames; //file name storage
  char* afname = (char*)malloc(bdim*sizeof(char)); //adapt file name storage
  char buff[bdim]; //buffer
  //char sname[bdim]; //mesh name storage
  time_t tm;
  char *t_char;
  double smn, smx, ar, ctol, mtol, nspace;
  int *boxcut = NULL;
  int *vbd = NULL;
  int *regions = NULL;
 
  //Start up MPI
  //MPI_Init(&argcs, &pArgs);
   
  //Find out process rank of current instance
  //MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
   
  //Find out number of processes
  //MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  
  //since not using MPI, set globals
  num_procs = 1;
  my_rank = 0;
  
  //create geometry (if not used as a library...normally already done)
  geometry* geom = new geometry();
  
  in_f = stdin;
  out_f = stdout;

  if (my_rank == 0)
  {
  
    //create journal file
    if ((jou_f=fopen("FASTAR.jou","w")) == NULL)
    {
      printf("\nCouldn't open FASTAR journal file");
      //MPI_Abort(MPI_COMM_WORLD,0);
      exit(0);
    }
  
    //check for standard input (if not used as a library)
    if (--argcs < 1)
    {
      printf("\nNo input file specified!");
      printf("\nUsing standard input!");
    } else
    {
      if ((in_f=fopen(pArgs[argcs],"r")) == NULL)
      {
        fprintf(stderr,"\nCouldn't open file <%s>\n",pArgs[argcs]);
        fprintf(stderr,"\n\nUsage: FASTAR.XXX [batch_input_file]\n");
        //MPI_Abort(MPI_COMM_WORLD,0);
        exit(0);
      }
    }
  }

  if (my_rank == 0 && in_f != stdin)
  {
    sprintf(buff,"FASTAR.out");

    if ((out_f=fopen(buff,"w")) == NULL)
    {
      fprintf(stderr,"\nCouldn't open file output file %s",buff);
      if (my_rank == 0)
        fclose(jou_f);
      //MPI_Abort(MPI_COMM_WORLD,0);
      exit(0);
    }
  }
  
  #ifdef _DEBUG
  buff[0] = '\0';
  
  sprintf(buff,"debug.%d",my_rank);

  if ((debug_f=fopen(buff,"w")) == NULL)
  {
    fprintf(stderr,"\nCouldn't open debug output file %s",buff);
    fflush(stderr);
    fclose(jou_f);
    fclose(out_f);
    //MPI_Abort(MPI_COMM_WORLD,0);
    exit(0);
  }
  #endif

  //MPI_Barrier(MPI_COMM_WORLD);

  //perform time check (if not used as a library)
  time(&tm);
  t_char = ctime(&tm);
  if (my_rank == 0)
  {
    fprintf(out_f,"\nFASTAR run started at %s",t_char);

    //print program info 
    fprintf(out_f,"\n===========================================================================");
    fprintf(out_f,"\n COPYRIGHT 2003-2010 THE UNIVERSITY OF TENNESSEE AT CHATTANOOGA            ");
    fprintf(out_f,"\n                                                                           ");
    fprintf(out_f,"\n                    RIGHTS IN DATA                                         ");
    fprintf(out_f,"\n                                                                           ");
    fprintf(out_f,"\n THIS SOFTWARE IS SUBMITTED WITH RESTRICTED RIGHTS UNDER GOVERNMENT        ");
    fprintf(out_f,"\n    CONTRACTS. USE, REPRODUCTION, OR DISCLOSURE IS SUBJECT TO              ");
    fprintf(out_f,"\n         RESTRICTIONS SET FORTH IN THESE CONTRACTS AND FEDERAL             ");
    fprintf(out_f,"\n              RIGHTS IN DATA CONTRACT CLAUSES.                             ");
    fprintf(out_f,"\n       ALL RIGHTS NOT RESERVED FOR THE GOVERNMENT ARE RETAINED BY          ");
    fprintf(out_f,"\n              THE UNIVERSITY OF TENNESSEE AT CHATTANOOGA                   ");
    fprintf(out_f,"\n                                                                           ");
    fprintf(out_f,"\n Function Adaptive Split Tree Anisotropic Refinement (FASTAR)              ");
    fprintf(out_f,"\n NOTE: This data includes the UT SimCenter at Chattanooga P_OPT, P_VLI, and");
    fprintf(out_f,"\n P_REFINE codes, which were developed under private non-government funding.");
    fprintf(out_f,"\n This software is submitted with limited rights to use, reproduce,         ");
    fprintf(out_f,"\n and disclose this data for Government Purposes only.                      ");
    fprintf(out_f,"\n Requests for access to the software for non-governmental purposes         ");
    fprintf(out_f,"\n should be referrred to                                                    "); 
    fprintf(out_f,"\n                                                                           ");
    fprintf(out_f,"\n    Dr. Steve Karman                  Mr. Vincent Betro                    "); 
    fprintf(out_f,"\n    Steve-Karman@utc.edu              Vincent-Betro@utc.edu                "); 
    fprintf(out_f,"\n    423-425-5492  or  423-425-5470    423-425-5434                         "); 
    fprintf(out_f,"\n                                                                           ");
    fprintf(out_f,"\n    University of Tennessee SimCenter at Chattanooga                       "); 
    fprintf(out_f,"\n    701 East M. L. King Boulevard                                          "); 
    fprintf(out_f,"\n    Chattanooga, TN 37403                                                  "); 
    fprintf(out_f,"\n===========================================================================\n");
    fflush(out_f);
  }

  //MPI_Barrier(MPI_COMM_WORLD);

  if (my_rank == 0)
  {
    //read number of geom files (if not used as a library...normally already done)
    journal(in_f, out_f, jou_f, "#Number of geometry files ->",nf);
    //check for adequate number of geom files
    if (nf <= 0)
    {
      fprintf(out_f,"\nNumber of geometry files must be > 0\n");
      fflush(out_f);
      //MPI_Abort(MPI_COMM_WORLD,0);
      exit(0);
    }
	  
    // inform slave processes
    //tag = 0;
    //for (dest = 1; dest < num_procs; dest++)
      //MPI_Send(&nf, 1, MPI_INT, dest, tag, MPI_COMM_WORLD);

    //read geometry (if not used as a library...normally already done)
    fnames = (char**)malloc(nf*sizeof(char*));
  
    for (i=0; i < nf; i++)
    {
      fnames[i] = (char*)malloc(bdim*sizeof(char));
      sprintf(buff,"#File %d ->",i+1);
      journal(in_f, out_f, jou_f, buff, fnames[i]);
      // inform slave processes of file names, using strlen for count
      //tag = i+1;
      //for (dest = 1; dest < num_procs; dest++)
        //MPI_Send(fnames[i], strlen(fnames[i])+1, MPI_CHAR, dest, tag, MPI_COMM_WORLD);
    }
    
    //read adapt files
    journal(in_f, out_f, jou_f, "#Enter 0 if no adaptation spacing file present, 1 if present ->",adapt);
    	  
    // inform slave processes
    //tag = 0;
    //for (dest = 1; dest < num_procs; dest++)
      //MPI_Send(&adapt, 1, MPI_INT, dest, tag, MPI_COMM_WORLD);
    if (adapt)
      {
      //read file name
      sprintf(buff,"#File ->");
      journal(in_f, out_f, jou_f, buff, afname);
      // inform slave processes of file name, using strlen for count
      //tag = 1;
      //for (dest = 1; dest < num_procs; dest++)
        //MPI_Send(afname, strlen(afname)+1, MPI_CHAR, dest, tag, MPI_COMM_WORLD); 
      }
  } else
  {
    // Obtain number of geometry files from master
    //source = 0; //reset to master
    //tag = 0;
    //MPI_Recv(&nf, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status);

    //read geometry (if not used as a library...normally already done)
    fnames = (char**)malloc(nf*sizeof(char*));
  
    for (i=0; i < nf; i++)
    {
      fnames[i] = (char*)malloc(bdim*sizeof(char));

      //tag = i+1; //based on master
      //MPI_Recv(fnames[i], bdim, MPI_CHAR, source, tag, MPI_COMM_WORLD, &status);
    }
    
    // Obtain adaptation file from master
    //source = 0; //reset to master
    //tag = 0;
    //MPI_Recv(&adapt, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
    
    if (adapt)
      {
      //read adapt
      //tag = 1; //based on master
      //MPI_Recv(afname, bdim, MPI_CHAR, source, tag, MPI_COMM_WORLD, &status);
      }
  }
  
  geom->read_geom(nf, fnames);
  
  //set up lo hi from geom for SF
  if (adapt)
    {
    Point l1, h1;
  
    l1 = geom->min_point();
    h1 = geom->max_point();
  
    double lo[3], hi[3];
  
    for (i = 0; i < 3; i++)
      {
      lo[i] = l1[i];
      hi[i] = h1[i];
      }

    SF_Initialize(afname, 1, 1, lo, hi);
    }
  
  //clean up fnames, afname
  for (i=0; i < nf; i++)
    free(fnames[i]);
  free(fnames);
  free(afname); //alloc'd regardless of adapt's value
  
  if (my_rank == 0)
  {  
    journal(in_f, out_f, jou_f,"#Enter 0 if geometry is to be used in refinement, 1 if only spacing file to be considered ->", ugeom);
    journal(in_f, out_f, jou_f,"#Enter minimum spacing ->", smn);
    journal(in_f, out_f, jou_f,"#Enter maximum spacing ->", smx);
    journal(in_f, out_f, jou_f,"#Enter maximum aspect ratio ->", ar);
    journal(in_f, out_f, jou_f,"#Enter 0 for speed, 1 for optimum quality ->", highquality);
    journal(in_f, out_f, jou_f,"#Enter 0 for standard flood-fill, 1 for nabor table ->", altflood);
    journal(in_f, out_f, jou_f,"#Enter number of test ray traces (will be cubed) ->", floodnum);
    journal(in_f, out_f, jou_f,"#Enter minimum metric length at which to refine (>=1.0) ->", mtol);
    journal(in_f, out_f, jou_f,"#Enter minimum spacing between cut surface and geometry ->", ctol);
    journal(in_f, out_f, jou_f,"#Enter desired normal spacing (be sure layers set to 0) ->", nspace);
    journal(in_f, out_f, jou_f,"#Enter 1 for restart, 0 for initial generation ->", restart);
    journal(in_f, out_f, jou_f,"#Enter 0 for UCD, 1 for STL ->", ftype);
    journal(in_f, out_f, jou_f,"#Enter 0 for normal cutting, 1 for box cutting ->", evencut);
    journal(in_f, out_f, jou_f,"#Enter number of boundaries to be box cut ->", bb);
    //allocate to hold bds
    boxcut = new int[bb];
    for (i = 0; i < bb; i++)
      journal(in_f, out_f, jou_f,"#Enter boundary number to be box cut ->", boxcut[i]);
    journal(in_f, out_f, jou_f,"#Enter 0 for output, 1 for TETGEN ->", gentet);
    journal(in_f, out_f, jou_f,"#Enter number of regions to be saved from Tetgen ->", numregions);
    regions = new int[numregions];
    for (i = 0; i < numregions; i++)
      journal(in_f, out_f, jou_f,"#Enter number of region to be stitched ->", regions[i]);
    journal(in_f, out_f, jou_f,"#Enter 0 for normal mode, 1 for overset ->", overset);
    journal(in_f, out_f, jou_f,"#Enter 0 for inviscid or a number of layers for viscous ->", viscous);
    if (viscous)
      {
      journal(in_f, out_f, jou_f,"#Enter geometric progression factor ->", factor);
      journal(in_f, out_f, jou_f,"#Enter intial viscous spacing ->", vspace);
      }
    journal(in_f, out_f, jou_f,"#Enter number of boundaries to be marked viscous ->", nvbd);
    //allocate to hold vbds
    vbd = new int[nvbd];
    for (i = 0; i < nvbd; i++)
      journal(in_f, out_f, jou_f,"#Enter boundary number to be marked viscous ->", vbd[i]);

    if (!ugeom)
      fprintf(out_f,"\nGeometry spacing file to be considered.");
    else
      fprintf(out_f,"\nOnly spacing/adaptation file to be considered.");
    fprintf(out_f,"\nMinimum spacing = %lg",smn);
    fprintf(out_f,"\nMaximum spacing = %lg",smx);
    fprintf(out_f,"\nMaximum aspect ratio = %lg",ar);
    if (!highquality)
      fprintf(out_f,"\nGenerating for speed.");
    else
      fprintf(out_f,"\nGenerating for optimum quality.");
    fprintf(out_f,"\nMinimum metric length at which to refine (>=1.0) = %lg",mtol);
    fprintf(out_f,"\nMinimum spacing between cut surface and geometry = %lg",ctol);
    fprintf(out_f,"\nDesired normal spacing (be sure layers set to 0) = %lg",nspace);
    fprintf(out_f,"\nRestart flag = %d",restart);
    if (!ftype)
      fprintf(out_f,"\nWill output UCD file for tet meshing.");
    else
      fprintf(out_f,"\nWill output STL file for tet meshing.");
    if (!evencut)
      fprintf(out_f,"\nWill do normal cutting.");
    else
      {
      fprintf(out_f,"\nWill do box cutting.");
      fprintf(out_f,"\nNumber of boundaries to box cut = %d.",bb);
      for (i = 0; i < bb; i++)
        fprintf(out_f,"\nBoundary %d will be box cut.",boxcut[i]);
      }
    if (gentet)
      fprintf(out_f,"\nWill use TetGen!");
    if (numregions)
      fprintf(out_f,"\nWill use %d regions from Tetgen.",numregions);
    for (i = 0; i < numregions; i++)
      fprintf(out_f,"\nWill use region %d from Tetgen.",regions[i]);
    if (altflood)
      fprintf(out_f,"\nWill do nabor table flood fill!");
    if (floodnum)
      fprintf(out_f,"\nWill do %d ray traces before flood fill!",floodnum*floodnum*floodnum);
    if (overset)
      fprintf(out_f,"\nWill create overset mesh.",overset);
    if (viscous)
      {
      fprintf(out_f,"\nWill insert %d viscous layers.",viscous);
      fprintf(out_f,"\nGeometric progression factor = %lf.",factor);
      fprintf(out_f,"\nInitial viscous spacing = %lf.",vspace);
      fprintf(out_f,"\nNumber of boundaries to be marked viscous = %d.",nvbd);
      for (i = 0; i < nvbd; i++)
        fprintf(out_f,"\nBoundary %d will be box cut.",vbd[i]);
      }
    fflush(out_f);
    
  
    //close journal file and input deck (if not used as a library)
    fclose(jou_f);
    if (in_f != stdin) fclose(in_f);
  }
  
  //inform all of min and max spacing
  //MPI_Bcast(&smn, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  //MPI_Bcast(&smx, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  //MPI_Bcast(&ar, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  //MPI_Bcast(&ctol, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  //MPI_Bcast(&mtol, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  //MPI_Bcast(&nspace, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  //MPI_Bcast(&restart, 1, MPI_INT, 0, MPI_COMM_WORLD);
  //MPI_Bcast(&ugeom, 1, MPI_INT, 0, MPI_COMM_WORLD);
  //MPI_Bcast(&ftype, 1, MPI_INT, 0, MPI_COMM_WORLD);
  //MPI_Bcast(&evencut, 1, MPI_INT, 0, MPI_COMM_WORLD);
  //MPI_Bcast(&altflood, 1, MPI_INT, 0, MPI_COMM_WORLD);
  //MPI_Bcast(&floodnum, 1, MPI_INT, 0, MPI_COMM_WORLD);
  //MPI_Bcast(&bb, 1, MPI_INT, 0, MPI_COMM_WORLD);
  //if (my_rank != 0)
    //boxcut = new int[bb];
  //MPI_Bcast(boxcut,bb,MPI_INT,0,MPI_COMM_WORLD);
  //MPI_Bcast(&gentet, 1, MPI_INT, 0, MPI_COMM_WORLD);
  //MPI_Bcast(&viscous, 1, MPI_INT, 0, MPI_COMM_WORLD);
  //MPI_Bcast(&factor, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  //MPI_Bcast(&vspace, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  //MPI_Bcast(&nvbd, 1, MPI_INT, 0, MPI_COMM_WORLD);
  //if (my_rank != 0)
    //vbd = new int[nvbd];
  //MPI_Bcast(vbd,nvbd,MPI_INT,0,MPI_COMM_WORLD);
  //MPI_Bcast(numregions,1,MPI_INT,0,MPI_COMM_WORLD);
  //if (my_rank != 0)
    //regions = new int[numregions];
  //MPI_Bcast(regions,numregions,MPI_INT,0,MPI_COMM_WORLD);

  //run p_opt as if it were a library with input parameters
  fastar_lib(geom, smn, smx, ar, ctol, mtol, nspace, boxcut, vbd, regions);
 
  //perform end time check
  time(&tm);
  t_char = ctime(&tm);
  if (my_rank == 0)
    {
    fprintf(out_f,"\nFASTAR run completed at %s",t_char);
    fflush(out_f);
    }
  
  //MPI_Barrier(MPI_COMM_WORLD);

  delete geom;

  if (boxcut != 0)
    delete [] boxcut;

  if (vbd != 0)
    delete [] vbd;

  if (regions != 0)
    delete [] regions;

  if (my_rank == 0 && out_f != stdout)
    fclose(out_f);
    
  #ifdef _DEBUG
  if (in_f != stdin)
    fclose(debug_f);
  #endif

  //MPI_Finalize(); 
  
  return(0);
}

