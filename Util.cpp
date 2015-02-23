#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Point.h"
#include "Vector.h"
#include "Util.h"

#define ABS(x) ((x) >= 0 ? (x) : -(x))
#define SIGN(a,b) ((b) < 0.0 ? -ABS(a) : ABS(a))

//global variables
extern int my_rank; /*rank of process*/
extern int num_procs; /*number of processes*/
extern FILE *in_f, *jou_f, *out_f, *debug_f; /*input output journal files global*/

int invert3X3(double m[3][3])
{
  double a, b, c, d, e, f, g, h, i;
  int flag;

  // invert matrix using adjoint via cofactors
  a = m[0][0];
  b = m[0][1];
  c = m[0][2];
  d = m[1][0];
  e = m[1][1];
  f = m[1][2];
  g = m[2][0];
  h = m[2][1];
  i = m[2][2];
  double det = a*e*i+b*f*g+c*d*h-c*e*g-b*d*i-a*f*h;
  if (fabs(det) > 1.0e-15)
    flag=1;
  else
    flag=0;
  det = SIGN(MAX(fabs(det),1.0e-15),det);
  m[0][0] =  (e*i-f*h)/det;
  m[0][1] = -(b*i-c*h)/det;
  m[0][2] =  (b*f-c*e)/det;
  m[1][0] = -(d*i-f*g)/det;
  m[1][1] =  (a*i-c*g)/det;
  m[1][2] = -(a*f-c*d)/det;
  m[2][0] =  (d*h-e*g)/det;
  m[2][1] = -(a*h-b*g)/det;
  m[2][2] =  (a*e-b*d)/det;
  
  return(flag);
}

double determinant_4X4(double m[4][4])
{
  double det;
  double m1[3][3], m2[3][3], m3[3][3], m4[3][3];

  m1[0][0] = m[1][1];
  m1[0][1] = m[1][2];
  m1[0][2] = m[1][3];
  m1[1][0] = m[2][1];
  m1[1][1] = m[2][2];
  m1[1][2] = m[2][3];
  m1[2][0] = m[3][1];
  m1[2][1] = m[3][2];
  m1[2][2] = m[3][3];
  m2[0][0] = m[0][1];
  m2[0][1] = m[0][2];
  m2[0][2] = m[0][3];
  m2[1][0] = m[2][1];
  m2[1][1] = m[2][2];
  m2[1][2] = m[2][3];
  m2[2][0] = m[3][1];
  m2[2][1] = m[3][2];
  m2[2][2] = m[3][3];
  m3[0][0] = m[0][1];
  m3[0][1] = m[0][2];
  m3[0][2] = m[0][3];
  m3[1][0] = m[1][1];
  m3[1][1] = m[1][2];
  m3[1][2] = m[1][3];
  m3[2][0] = m[3][1];
  m3[2][1] = m[3][2];
  m3[2][2] = m[3][3];
  m4[0][0] = m[0][1];
  m4[0][1] = m[0][2];
  m4[0][2] = m[0][3];
  m4[1][0] = m[1][1];
  m4[1][1] = m[1][2];
  m4[1][2] = m[1][3];
  m4[2][0] = m[2][1];
  m4[2][1] = m[2][2];
  m4[2][2] = m[2][3];

  det = m[0][0]*(m1[0][0]*m1[1][1]*m1[2][2]+
                 m1[0][1]*m1[1][2]*m1[2][0]+
                 m1[0][2]*m1[1][0]*m1[2][1]-
                 m1[0][2]*m1[1][1]*m1[2][0]-
                 m1[0][1]*m1[1][0]*m1[2][2]-
                 m1[0][0]*m1[1][2]*m1[2][1])-
        m[1][0]*(m2[0][0]*m2[1][1]*m2[2][2]+
                 m2[0][1]*m2[1][2]*m2[2][0]+
                 m2[0][2]*m2[1][0]*m2[2][1]-
                 m2[0][2]*m2[1][1]*m2[2][0]-
                 m2[0][1]*m2[1][0]*m2[2][2]-
                 m2[0][0]*m2[1][2]*m2[2][1])+
        m[2][0]*(m3[0][0]*m3[1][1]*m3[2][2]+
                 m3[0][1]*m3[1][2]*m3[2][0]+
                 m3[0][2]*m3[1][0]*m3[2][1]-
                 m3[0][2]*m3[1][1]*m3[2][0]-
                 m3[0][1]*m3[1][0]*m3[2][2]-
                 m3[0][0]*m3[1][2]*m3[2][1])-
        m[3][0]*(m4[0][0]*m4[1][1]*m4[2][2]+
                 m4[0][1]*m4[1][2]*m4[2][0]+
                 m4[0][2]*m4[1][0]*m4[2][1]-
                 m4[0][2]*m4[1][1]*m4[2][0]-
                 m4[0][1]*m4[1][0]*m4[2][2]-
                 m4[0][0]*m4[1][2]*m4[2][1]);

  return(det);
}

double distance( Point p1, Point p2)
{
  double dx = p2[0]-p1[0];
  double dy = p2[1]-p1[1];
  double dz = p2[2]-p1[2];
  double mag = sqrt(dx*dx+dy*dy+dz*dz);
  return mag;
}

double closest_to_edge(Point p1, Point p2, Point p, Point &ep)
{
  double mag, d1, d2, dot;
  Vector v1, v2;

  v1 = Vector(p1,p2);
  v2 = Vector(p1,p);
  d1 = v1*v1;
  dot = (v1*v2)/MAX(1.0e-15,d1);
  d2 = MIN(1.0,MAX(0.0,dot));
  d1 = 1.0-d2;
  ep = p1*d1+p2*d2;
  v1 = Vector(p,ep);
  mag = v1.magnitude();
  return(mag);
}

double closest_to_facet(Point p0, Point p1, Point p2, Point p, Point &pt)
{
  double mag, d, dot, a0, a1, a2;
  Vector v1, v2, v3, norm, area;
  double tol = 1.0e-15;

  v1 = Vector(p0,p1);
  v2 = Vector(p0,p2);
  norm = (v1 % v2);
  mag = MAX(tol,norm.magnitude());
  //norm.normalize();
  norm /= mag;
  v1 = Vector(p0,p);
  d = v1*norm;
  v1 -= norm*d;
  pt[0] = p0[0] + v1[0];
  pt[1] = p0[1] + v1[1];
  pt[2] = p0[2] + v1[2];
  v1 = Vector(pt,p0);
  v2 = Vector(pt,p1);
  v3 = Vector(pt,p2);
  area = (v2 % v3);
  a0 = (area*norm)/mag;
  area = (v3 % v1);
  a1 = (area*norm)/mag;
  area = (v1 % v2);
  a2 = (area*norm)/mag;
  // Restrict test point to within triangle boundaries
  if (a0 < 0.0)
  {
    v1 = Vector(p1,p);
    v2 = Vector(p1,p2);
    dot = (v1*v2)/MAX(tol,v2*v2);
    a2 = MIN(1.0,MAX(0.0,dot));
    a1 = 1.0-a2;
    pt = p1*a1+p2*a2;
  } else if (a1 < 0.0)
  {
    v1 = Vector(p2,p);
    v2 = Vector(p2,p0);
    dot = (v1*v2)/MAX(tol,v2*v2);
    a0 = MIN(1.0,MAX(0.0,dot));
    a2 = 1.0-a0;
    pt = p0*a0+p2*a2;
  } else if (a2 < 0.0)
  {
    v1 = Vector(p0,p);
    v2 = Vector(p0,p1);
    dot = (v1*v2)/MAX(tol,v2*v2);
    a1 = MIN(1.0,MAX(0.0,dot));
    a0 = 1.0-a1;
    pt = p0*a0+p1*a1;
  }
  mag = distance(pt,p);

  return(mag);
}

int line_line_intersect(Point pa, Point pb, Point p1, Point p2, Point &p)
{

  double rdnm, tol, mag, dab, d12;
 
  double eps = 1.0e-5;
  Vector v1, v2;
  v1 = Vector(pa,pb);
  v2 = Vector(p1,p2);
  tol = MIN(v1.magnitude(),v2.magnitude())*1.0e-15;
  tol = MAX(tol,1.0e-15);
  rdnm = (pa[1]-pb[1])*(p2[0]-p1[0])-(pb[0]-pa[0])*(p1[1]-p2[1]);
  if (fabs(rdnm) < tol) 
    return(0);
  rdnm = 1.0/rdnm;
  p[0] = ((pb[0]*pa[1]-pa[0]*pb[1])*(p2[0]-p1[0])-
          (pb[0]-pa[0])*(p2[0]*p1[1]-p1[0]*p2[1]))*rdnm;
  p[1] = ((pa[1]-pb[1])*(p2[0]*p1[1]-p1[0]*p2[1])-
          (pb[0]*pa[1]-pa[0]*pb[1])*(p1[1]-p2[1]))*rdnm;
  v1 = Vector(pa,pb);
  v2 = Vector(pa,p);
  mag = v1 * v1;
  dab = (v1 * v2)/MAX(tol,mag);
  v1 = Vector(p1,p2);
  v2 = Vector(p1,p);
  mag = v1 * v1;
  d12 = (v1 * v2)/MAX(tol,mag);
  if ((dab >= -eps && dab <= 1.0+eps) && (d12 >= -eps && d12 <= 1.0+eps))
    return(1);
  else
    return(0);
}

int facet_facet_intersect(Point p1, Point p2, Point p3, Point q1, Point q2, Point q3, Point &pa, Point &pb, double tol)
{
  Vector edge, line, pnorm, qnorm, v1, v2;
  double dp, dq, emag, pmag, qmag, rdnm, s;
  Point lpt, pi[6], qi[6], t0, t1, ps;
  int i, j, npi, nqi;
  double ptol = 1.0-1.0e-12;
  double dtol = 1.0e-15;
  double etol = 0.0e-05;
  FILE *fp;
  char buff[132];

  // normal vector & plane equation for plane "p"
  v1 = Vector(p1, p2);
  v2 = Vector(p1, p3);
  pnorm = v1 % v2;
  pmag = MAX(dtol,pnorm.magnitude());
  dp = -(p1[0]*p2[1]*p3[2] + p1[1]*p2[2]*p3[0] + p1[2]*p2[0]*p3[1]
        -p1[2]*p2[1]*p3[0] - p1[1]*p2[0]*p3[2] - p1[0]*p2[2]*p3[1]);
  pnorm /= pmag;
  dp /= pmag;

  // normal vector & plane equation for plane "q"
  v1 = Vector(q1, q2);
  v2 = Vector(q1, q3);
  qnorm = v1 % v2;
  qmag = MAX(dtol,qnorm.magnitude());
  dq = -(q1[0]*q2[1]*q3[2] + q1[1]*q2[2]*q3[0] + q1[2]*q2[0]*q3[1]
        -q1[2]*q2[1]*q3[0] - q1[1]*q2[0]*q3[2] - q1[0]*q2[2]*q3[1]);
  qnorm /= qmag;
  dq /= qmag;

  if (fabs(pnorm * qnorm) > ptol)
    return(0);

  // direction of intersection line is perpendicular to p & q
  line = pnorm % qnorm;
  line.normalize();

  // find a point on the intersection line to use for parametric coordinate calculation
  if (fabs(line[2]) >= MAX(fabs(line[0]),fabs(line[1])))
  {
    lpt[2] = 0.0;
    rdnm = pnorm[0]*qnorm[1]-pnorm[1]*qnorm[0];
    rdnm = 1.0/SIGN(MAX(1.0e-15,fabs(rdnm)),rdnm);
    lpt[0] = ((-dp)*qnorm[1] - pnorm[1]*(-dq))*rdnm;
    lpt[1] = (pnorm[0]*(-dq) - (-dp)*qnorm[0])*rdnm;
  } else if (fabs(line[1]) >= MAX(fabs(line[0]),fabs(line[2])))
  {
    lpt[1] = 0.0;
    rdnm = pnorm[0]*qnorm[2]-pnorm[2]*qnorm[0];
    rdnm = 1.0/SIGN(MAX(1.0e-15,fabs(rdnm)),rdnm);
    lpt[0] = ((-dp)*qnorm[2] - pnorm[2]*(-dq))*rdnm;
    lpt[2] = (pnorm[0]*(-dq) - (-dp)*qnorm[0])*rdnm;
  } else if (fabs(line[0]) >= MAX(fabs(line[1]),fabs(line[2])))
  {
    lpt[0] = 0.0;
    rdnm = pnorm[1]*qnorm[2]-pnorm[2]*qnorm[1];
    rdnm = 1.0/SIGN(MAX(1.0e-15,fabs(rdnm)),rdnm);
    lpt[1] = ((-dp)*qnorm[2] - pnorm[2]*(-dq))*rdnm;
    lpt[2] = (pnorm[1]*(-dq) - (-dp)*qnorm[1])*rdnm;
  } else
  {
    fprintf(stderr,"facet_facet_intersection: Unable to find point on line!");
    exit(0);
  }

  // now test intersection of define line and edges of each facet

  npi = nqi = 0;

  for (i=0; i < 3; i++)
  {
    switch(i)
    {
      case 0:
        t0 = p1;
        t1 = p2;
        break;
      case 1:
        t0 = p2;
        t1 = p3;
        break;
      case 2:
        t0 = p3;
        t1 = p1;
        break;
    }
    edge = Vector(t0,t1);
    emag = edge.magnitude();
    edge /= MAX(dtol,emag);
    if (fabs(edge * line) < ptol) // line and edge not parallel
    {
      // choose the Cartesian plane of the largest normal component
      if (fabs(pnorm[2]) >= MAX(fabs(pnorm[0]),fabs(pnorm[1])))
      {
        rdnm = -edge[0]*line[1] + line[0]*edge[1];
        rdnm = 1.0/SIGN(MAX(dtol,fabs(rdnm)),rdnm);
        s = (-(lpt[0]-t0[0])*line[1] + line[0]*(lpt[1]-t0[1]))*rdnm;
        s /= emag;
        if (s >= -etol && s <= 1.0+etol)
          pi[npi++] = t0 + (t1-t0)*s;
      } else if (fabs(pnorm[1]) >= MAX(fabs(pnorm[0]),fabs(pnorm[2])))
      {
        rdnm = -edge[0]*line[2] + line[0]*edge[2];
        rdnm = 1.0/SIGN(MAX(dtol,fabs(rdnm)),rdnm);
        s = (-(lpt[0]-t0[0])*line[2] + line[0]*(lpt[2]-t0[2]))*rdnm;
        s /= emag;
        if (s >= -etol && s <= 1.0+etol)
          pi[npi++] = t0 + (t1-t0)*s;
      } else
      {
        rdnm = -edge[1]*line[2] + line[1]*edge[2];
        rdnm = 1.0/SIGN(MAX(dtol,fabs(rdnm)),rdnm);
        s = (-(lpt[1]-t0[1])*line[2] + line[1]*(lpt[2]-t0[2]))*rdnm;
        s /= emag;
        if (s >= -etol && s <= 1.0+etol)
          pi[npi++] = t0 + (t1-t0)*s;
      }
    } else
    {
      // lines parallel, check for co-linear
      s = line * Vector(lpt,t0);
      ps= lpt + Point(line[0],line[1],line[2])*s;
      if (distance(ps,t0) <= tol)
      {
        pi[npi++] = t0;
        pi[npi++] = t1;
      }
    }

    switch(i)
    {
      case 0:
        t0 = q1;
        t1 = q2;
        break;
      case 1:
        t0 = q2;
        t1 = q3;
        break;
      case 2:
        t0 = q3;
        t1 = q1;
        break;
    }
    edge = Vector(t0,t1);
    emag = edge.magnitude();
    edge /= MAX(dtol,emag);
    if (fabs(edge * line) < ptol) // line and edge not parallel
    {
      // choose the Cartesian plane of the largest normal component
      if (fabs(qnorm[2]) >= MAX(fabs(qnorm[0]),fabs(qnorm[1])))
      {
        rdnm = -edge[0]*line[1] + line[0]*edge[1];
        rdnm = 1.0/SIGN(MAX(dtol,fabs(rdnm)),rdnm);
        s = (-(lpt[0]-t0[0])*line[1] + line[0]*(lpt[1]-t0[1]))*rdnm;
        s /= emag;
        if (s >= -etol && s <= 1.0+etol)
          qi[nqi++] = t0 + (t1-t0)*s;
      } else if (fabs(qnorm[1]) >= MAX(fabs(qnorm[0]),fabs(qnorm[2])))
      {
        rdnm = -edge[0]*line[2] + line[0]*edge[2];
        rdnm = 1.0/SIGN(MAX(dtol,fabs(rdnm)),rdnm);
        s = (-(lpt[0]-t0[0])*line[2] + line[0]*(lpt[2]-t0[2]))*rdnm;
        s /= emag;
        if (s >= -etol && s <= 1.0+etol)
          qi[nqi++] = t0 + (t1-t0)*s;
      } else
      {
        rdnm = -edge[1]*line[2] + line[1]*edge[2];
        rdnm = 1.0/SIGN(MAX(dtol,fabs(rdnm)),rdnm);
        s = (-(lpt[1]-t0[1])*line[2] + line[1]*(lpt[2]-t0[2]))*rdnm;
        s /= emag;
        if (s >= -etol && s <= 1.0+etol)
          qi[nqi++] = t0 + (t1-t0)*s;
      }
    } else
    {
      // lines parallel, check for co-linear
      s = line * Vector(lpt,t0);
      ps= lpt + Point(line[0],line[1],line[2])*s;
      if (distance(ps,t0) <= tol)
      {
        qi[nqi++] = t0;
        qi[nqi++] = t1;
      }
    }
  }

  double ds, dsmn, ltol;

  // check for duplicate saved points
  ltol = tol;
  do
  {
    dsmn = 1.0e20;
    for (i=0; i < npi-1; i++)
    {
      j = i+1;
      while (j < npi)
      {
        ds = distance(pi[i],pi[j]);
        dsmn = MIN(ds,dsmn);
        if (ds < ltol)
          pi[j] = pi[--npi];
        else
          j++;
      }
    }
    ltol = dsmn*1.01;
  } while (npi > 2);

  ltol = tol;
  do
  {
    for (i=0; i < nqi-1; i++)
    {
      j = i+1;
      while (j < nqi)
      {
        ds = distance(qi[i],qi[j]);
        dsmn = MIN(ds,dsmn);
        if (ds < ltol)
          qi[j] = qi[--nqi];
        else
          j++;
      }
    }
    ltol = dsmn*1.01;
  } while (nqi > 2);

  if (npi <= 1 || nqi <= 1) // not enough intersection points on one or both facets
    return(0);

  if (npi > 2)
  {
    fprintf(stderr,"facet_facet_intersect: How can I have %d intersection points for plane p?",npi);
    sprintf(buff,"facet_p.dat");
    fp=fopen(buff,"w");
    fprintf(fp,"%14.7e %14.7e %14.7e\n",p1[0],p1[1],p1[2]);
    fprintf(fp,"%14.7e %14.7e %14.7e\n",p2[0],p2[1],p2[2]);
    fprintf(fp,"%14.7e %14.7e %14.7e\n",p3[0],p3[1],p3[2]);
    fprintf(fp,"%14.7e %14.7e %14.7e\n",p1[0],p1[1],p1[2]);
    fclose(fp);
    sprintf(buff,"facet_q.dat");
    fp=fopen(buff,"w");
    fprintf(fp,"%14.7e %14.7e %14.7e\n",q1[0],q1[1],q1[2]);
    fprintf(fp,"%14.7e %14.7e %14.7e\n",q2[0],q2[1],q2[2]);
    fprintf(fp,"%14.7e %14.7e %14.7e\n",q3[0],q3[1],q3[2]);
    fprintf(fp,"%14.7e %14.7e %14.7e\n",q1[0],q1[1],q1[2]);
    fclose(fp);
    sprintf(buff,"facet_pi.dat");
    fp=fopen(buff,"w");
    for (i=0; i < npi; i++)
      fprintf(fp,"%14.7e %14.7e %14.7e\n",pi[i][0],pi[i][1],pi[i][2]);
    fclose(fp);
    sprintf(buff,"facet_qi.dat");
    fp=fopen(buff,"w");
    for (i=0; i < nqi; i++)
      fprintf(fp,"%14.7e %14.7e %14.7e\n",qi[i][0],qi[i][1],qi[i][2]);
    fclose(fp);
    exit(0);
  }
  if (nqi > 2)
  {
    fprintf(stderr,"facet_facet_intersect: How can I have %d intersection points for plane q?",nqi);
    sprintf(buff,"facet_p.dat");
    fp=fopen(buff,"w");
    fprintf(fp,"%14.7e %14.7e %14.7e\n",p1[0],p1[1],p1[2]);
    fprintf(fp,"%14.7e %14.7e %14.7e\n",p2[0],p2[1],p2[2]);
    fprintf(fp,"%14.7e %14.7e %14.7e\n",p3[0],p3[1],p3[2]);
    fprintf(fp,"%14.7e %14.7e %14.7e\n",p1[0],p1[1],p1[2]);
    fclose(fp);
    sprintf(buff,"facet_q.dat");
    fp=fopen(buff,"w");
    fprintf(fp,"%14.7e %14.7e %14.7e\n",q1[0],q1[1],q1[2]);
    fprintf(fp,"%14.7e %14.7e %14.7e\n",q2[0],q2[1],q2[2]);
    fprintf(fp,"%14.7e %14.7e %14.7e\n",q3[0],q3[1],q3[2]);
    fprintf(fp,"%14.7e %14.7e %14.7e\n",q1[0],q1[1],q1[2]);
    fclose(fp);
    sprintf(buff,"facet_pi.dat");
    fp=fopen(buff,"w");
    for (i=0; i < npi; i++)
      fprintf(fp,"%14.7e %14.7e %14.7e\n",pi[i][0],pi[i][1],pi[i][2]);
    fclose(fp);
    sprintf(buff,"facet_qi.dat");
    fp=fopen(buff,"w");
    for (i=0; i < nqi; i++)
      fprintf(fp,"%14.7e %14.7e %14.7e\n",qi[i][0],qi[i][1],qi[i][2]);
    fclose(fp);
    exit(0);
  }

  // for single save nodes copy for sorting logic
  if (npi == 1)
    pi[npi++] = pi[0];
  if (nqi == 1)
    qi[nqi++] = qi[0];

  // sort the save points in each facet according to line direction
  Point swap;
  v1 = Vector(pi[0],pi[1]);
  if (v1 * line < 0.0)
  {
    swap = pi[0];
    pi[0] = pi[1];
    pi[1] = swap;
    v1 = Vector(pi[0],pi[1]);
  }
  v2 = Vector(qi[0],qi[1]);
  if (v2 * line < 0.0)
  {
    swap = qi[0];
    qi[0] = qi[1];
    qi[1] = swap;
    v2 = Vector(qi[0],qi[1]);
  }

  if ((line * Vector(lpt,pi[0])) > (line * Vector(lpt,qi[0])))
    pa = pi[0];
  else
    pa = qi[0];

  if ((line * Vector(lpt,pi[1])) < (line * Vector(lpt,qi[1])))
    pb = pi[1];
  else
    pb = qi[1];

  if ((line * Vector(lpt,pb)) < (line * Vector(lpt,pa)))
    return(0);
  else
  {
    // ensure segment is in direction of line vector
    v1 = Vector(pa,pb);
    if (v1 * line < 0.0)
    {
      swap = pa;
      pa = pb;
      pb = swap;
    }
#ifdef _DEBUG
//    sprintf(buff,"facet_segment.dat");
//    fp=fopen(buff,"w");
//    fprintf(fp,"%14.7e %14.7e %14.7e\n",pa[0],pa[1],pa[2]);
//    fprintf(fp,"%14.7e %14.7e %14.7e\n",pb[0],pb[1],pb[2]);
//    fclose(fp);
#endif
    return(1);
  }
  
}

int tri_quad_intersect(Point t1, Point t2, Point t3, Point q1, Point q2, Point q3, Point q4, Point &pa, Point &pb, double tol)
{
  Vector edge, line, tnorm, qnorm, v1, v2;
  double dt, dq, emag, tmag, qmag, rdnm, s;
  Point lpt, ti[6], qi[6], p0, p1, ps, lo, hi, mid;
  int i, j, nti, nqi;
  double ptol = 1.0-1.0e-12;
  double dtol = 1.0e-15;
  double etol = 1.0e-10;
  double stol = 1.0e-10;
  FILE *fp;
  char buff[132];

  // THIS ROUTINE ASSUMES THE QUAD IS PLANAR!

  lo = t1;
  hi = t1;
  for (i=0; i < 3; i++)
  {
    lo[i] = MIN(lo[i],t2[i]);
    lo[i] = MIN(lo[i],t3[i]);
    lo[i] = MIN(lo[i],q1[i]);
    lo[i] = MIN(lo[i],q2[i]);
    lo[i] = MIN(lo[i],q3[i]);
    lo[i] = MIN(lo[i],q4[i]);
    hi[i] = MAX(hi[i],t2[i]);
    hi[i] = MAX(hi[i],t3[i]);
    hi[i] = MAX(hi[i],q1[i]);
    hi[i] = MAX(hi[i],q2[i]);
    hi[i] = MAX(hi[i],q3[i]);
    hi[i] = MAX(hi[i],q4[i]);
  }
  mid = (lo + hi)*0.5;

  // normal vector & plane equation for triangle plane
  v1 = Vector(t1, t2);
  v2 = Vector(t1, t3);
  tnorm = v1 % v2;
  tmag = MAX(dtol,tnorm.magnitude());
  dt = -(t1[0]*t2[1]*t3[2] + t1[1]*t2[2]*t3[0] + t1[2]*t2[0]*t3[1]
        -t1[2]*t2[1]*t3[0] - t1[1]*t2[0]*t3[2] - t1[0]*t2[2]*t3[1]);
  tnorm /= tmag;
  dt /= tmag;

  // normal vector & plane equation for plane "q"
  v1 = Vector(q1, q2);
  v2 = Vector(q1, q3);
  qnorm = v1 % v2;
  qmag = MAX(dtol,qnorm.magnitude());
  dq = -(q1[0]*q2[1]*q3[2] + q1[1]*q2[2]*q3[0] + q1[2]*q2[0]*q3[1]
        -q1[2]*q2[1]*q3[0] - q1[1]*q2[0]*q3[2] - q1[0]*q2[2]*q3[1]);
  qnorm /= qmag;
  dq /= qmag;

  if (fabs(tnorm * qnorm) > ptol)
    return(0);

  // direction of intersection line is perpendicular to p & q
  line = tnorm % qnorm;
  line.normalize();

  // find a point on the intersection line to use for parametric coordinate calculation
  if (fabs(line[2]) >= MAX(fabs(line[0]),fabs(line[1])))
  {
    lpt[2] = mid[2];
    rdnm = tnorm[0]*qnorm[1]-tnorm[1]*qnorm[0];
    rdnm = 1.0/SIGN(MAX(1.0e-15,fabs(rdnm)),rdnm);
    lpt[0] = ((-dt-tnorm[2]*mid[2])*qnorm[1] - tnorm[1]*(-dq-qnorm[2]*mid[2]))*rdnm;
    lpt[1] = (tnorm[0]*(-dq-qnorm[2]*mid[2]) - (-dt-tnorm[2]*mid[2])*qnorm[0])*rdnm;
  } else if (fabs(line[1]) >= MAX(fabs(line[0]),fabs(line[2])))
  {
    lpt[1] = mid[1];
    rdnm = tnorm[0]*qnorm[2]-tnorm[2]*qnorm[0];
    rdnm = 1.0/SIGN(MAX(1.0e-15,fabs(rdnm)),rdnm);
    lpt[0] = ((-dt-tnorm[1]*mid[1])*qnorm[2] - tnorm[2]*(-dq-qnorm[1]*mid[1]))*rdnm;
    lpt[2] = (tnorm[0]*(-dq-qnorm[1]*mid[1]) - (-dt-tnorm[1]*mid[1])*qnorm[0])*rdnm;
  } else if (fabs(line[0]) >= MAX(fabs(line[1]),fabs(line[2])))
  {
    lpt[0] = mid[0];
    rdnm = tnorm[1]*qnorm[2]-tnorm[2]*qnorm[1];
    rdnm = 1.0/SIGN(MAX(1.0e-15,fabs(rdnm)),rdnm);
    lpt[1] = ((-dt-tnorm[0]*mid[0])*qnorm[2] - tnorm[2]*(-dq-qnorm[0]*mid[0]))*rdnm;
    lpt[2] = (tnorm[1]*(-dq-qnorm[0]*mid[0]) - (-dt-tnorm[0]*mid[0])*qnorm[1])*rdnm;
  } else
  {
    fprintf(stderr,"facet_facet_intersection: Unable to find point on line!");
    exit(0);
  }

  // now test intersection of define line and edges of each facet

  nti = nqi = 0;

  for (i=0; i < 3; i++)
  {
    switch(i)
    {
      case 0:
        p0 = t1;
        p1 = t2;
        break;
      case 1:
        p0 = t2;
        p1 = t3;
        break;
      case 2:
        p0 = t3;
        p1 = t1;
        break;
    }
    edge = Vector(p0,p1);
    emag = edge.magnitude();
    edge /= MAX(dtol,emag);
    if (fabs(edge * line) < ptol) // line and edge not parallel
    {
      // choose the Cartesian plane of the largest normal component
      if (fabs(tnorm[2]) >= MAX(fabs(tnorm[0]),fabs(tnorm[1])))
      {
        rdnm = -edge[0]*line[1] + line[0]*edge[1];
        rdnm = 1.0/SIGN(MAX(dtol,fabs(rdnm)),rdnm);
        s = (-(lpt[0]-p0[0])*line[1] + line[0]*(lpt[1]-p0[1]))*rdnm;
        s /= emag;
        if (s >= -etol && s <= 1.0+etol)
          ti[nti++] = p0 + (p1-p0)*s;
      } else if (fabs(tnorm[1]) >= MAX(fabs(tnorm[0]),fabs(tnorm[2])))
      {
        rdnm = -edge[0]*line[2] + line[0]*edge[2];
        rdnm = 1.0/SIGN(MAX(dtol,fabs(rdnm)),rdnm);
        s = (-(lpt[0]-p0[0])*line[2] + line[0]*(lpt[2]-p0[2]))*rdnm;
        s /= emag;
        if (s >= -etol && s <= 1.0+etol)
          ti[nti++] = p0 + (p1-p0)*s;
      } else
      {
        rdnm = -edge[1]*line[2] + line[1]*edge[2];
        rdnm = 1.0/SIGN(MAX(dtol,fabs(rdnm)),rdnm);
        s = (-(lpt[1]-p0[1])*line[2] + line[1]*(lpt[2]-p0[2]))*rdnm;
        s /= emag;
        if (s >= -etol && s <= 1.0+etol)
          ti[nti++] = p0 + (p1-p0)*s;
      }
    } else
    {
      // lines parallel, check for co-linear
      s = line * Vector(lpt,p0);
      ps= lpt + Point(line[0],line[1],line[2])*s;
      if (distance(ps,p0) <= stol)
      {
        ti[nti++] = p0;
        ti[nti++] = p1;
      }
    }
  }

  for (i=0; i < 4; i++)
  {
    switch(i)
    {
      case 0:
        p0 = q1;
        p1 = q2;
        break;
      case 1:
        p0 = q2;
        p1 = q3;
        break;
      case 2:
        p0 = q3;
        p1 = q4;
        break;
      case 3:
        p0 = q4;
        p1 = q1;
        break;
    }
    edge = Vector(p0,p1);
    emag = edge.magnitude();
    edge /= MAX(dtol,emag);
    if (fabs(edge * line) < ptol) // line and edge not parallel
    {
      // choose the Cartesian plane of the largest normal component
      if (fabs(qnorm[2]) >= MAX(fabs(qnorm[0]),fabs(qnorm[1])))
      {
        rdnm = -edge[0]*line[1] + line[0]*edge[1];
        rdnm = 1.0/SIGN(MAX(dtol,fabs(rdnm)),rdnm);
        s = (-(lpt[0]-p0[0])*line[1] + line[0]*(lpt[1]-p0[1]))*rdnm;
        s /= emag;
        if (s >= -etol && s <= 1.0+etol)
          qi[nqi++] = p0 + (p1-p0)*s;
      } else if (fabs(qnorm[1]) >= MAX(fabs(qnorm[0]),fabs(qnorm[2])))
      {
        rdnm = -edge[0]*line[2] + line[0]*edge[2];
        rdnm = 1.0/SIGN(MAX(dtol,fabs(rdnm)),rdnm);
        s = (-(lpt[0]-p0[0])*line[2] + line[0]*(lpt[2]-p0[2]))*rdnm;
        s /= emag;
        if (s >= -etol && s <= 1.0+etol)
          qi[nqi++] = p0 + (p1-p0)*s;
      } else
      {
        rdnm = -edge[1]*line[2] + line[1]*edge[2];
        rdnm = 1.0/SIGN(MAX(dtol,fabs(rdnm)),rdnm);
        s = (-(lpt[1]-p0[1])*line[2] + line[1]*(lpt[2]-p0[2]))*rdnm;
        s /= emag;
        if (s >= -etol && s <= 1.0+etol)
          qi[nqi++] = p0 + (p1-p0)*s;
      }
    } else
    {
      // lines parallel, check for co-linear
      s = line * Vector(lpt,p0);
      ps= lpt + Point(line[0],line[1],line[2])*s;
      if (distance(ps,p0) <= stol)
      {
        qi[nqi++] = p0;
        qi[nqi++] = p1;
      }
    }
  }

  double ds, dsmn, ltol;

  // check for duplicate saved points
  //ltol = tol;
  ltol = stol;
  do
  {
    dsmn = 1.0e20;
    for (i=0; i < nti-1; i++)
    {
      j = i+1;
      while (j < nti)
      {
        ds = distance(ti[i],ti[j]);
        dsmn = MIN(ds,dsmn);
        if (ds < ltol)
          ti[j] = ti[--nti];
        else
          j++;
      }
    }
    ltol = dsmn*1.01;
  } while (nti > 2);

  //ltol = tol;
  ltol = stol;
  do
  {
    for (i=0; i < nqi-1; i++)
    {
      j = i+1;
      while (j < nqi)
      {
        ds = distance(qi[i],qi[j]);
        dsmn = MIN(ds,dsmn);
        if (ds < ltol)
          qi[j] = qi[--nqi];
        else
          j++;
      }
    }
    ltol = dsmn*1.01;
  } while (nqi > 2);

  if (nti <= 1 || nqi <= 1) // not enough intersection points on one or both facets
    return(0);

  if (nti > 2)
  {
    fprintf(stderr,"tri_quad_intersect: How can I have %d intersection points for triangle plane?",nti);
    sprintf(buff,"facet_t.dat");
    fp=fopen(buff,"w");
    fprintf(fp,"%14.7e %14.7e %14.7e\n",t1[0],t1[1],t1[2]);
    fprintf(fp,"%14.7e %14.7e %14.7e\n",t2[0],t2[1],t2[2]);
    fprintf(fp,"%14.7e %14.7e %14.7e\n",t3[0],t3[1],t3[2]);
    fprintf(fp,"%14.7e %14.7e %14.7e\n",t1[0],t1[1],t1[2]);
    fclose(fp);
    sprintf(buff,"facet_q.dat");
    fp=fopen(buff,"w");
    fprintf(fp,"%14.7e %14.7e %14.7e\n",q1[0],q1[1],q1[2]);
    fprintf(fp,"%14.7e %14.7e %14.7e\n",q2[0],q2[1],q2[2]);
    fprintf(fp,"%14.7e %14.7e %14.7e\n",q3[0],q3[1],q3[2]);
    fprintf(fp,"%14.7e %14.7e %14.7e\n",q4[0],q4[1],q4[2]);
    fprintf(fp,"%14.7e %14.7e %14.7e\n",q1[0],q1[1],q1[2]);
    fclose(fp);
    sprintf(buff,"facet_ti.dat");
    fp=fopen(buff,"w");
    for (i=0; i < nti; i++)
      fprintf(fp,"%14.7e %14.7e %14.7e\n",ti[i][0],ti[i][1],ti[i][2]);
    fclose(fp);
    sprintf(buff,"facet_qi.dat");
    fp=fopen(buff,"w");
    for (i=0; i < nqi; i++)
      fprintf(fp,"%14.7e %14.7e %14.7e\n",qi[i][0],qi[i][1],qi[i][2]);
    fclose(fp);
    exit(0);
  }
  if (nqi > 2)
  {
    fprintf(stderr,"tri_quad_intersect: How can I have %d intersection points for quad plane?",nqi);
    sprintf(buff,"facet_t.dat");
    fp=fopen(buff,"w");
    fprintf(fp,"%14.7e %14.7e %14.7e\n",t1[0],t1[1],t1[2]);
    fprintf(fp,"%14.7e %14.7e %14.7e\n",t2[0],t2[1],t2[2]);
    fprintf(fp,"%14.7e %14.7e %14.7e\n",t3[0],t3[1],t3[2]);
    fprintf(fp,"%14.7e %14.7e %14.7e\n",t1[0],t1[1],t1[2]);
    fclose(fp);
    sprintf(buff,"facet_q.dat");
    fp=fopen(buff,"w");
    fprintf(fp,"%14.7e %14.7e %14.7e\n",q1[0],q1[1],q1[2]);
    fprintf(fp,"%14.7e %14.7e %14.7e\n",q2[0],q2[1],q2[2]);
    fprintf(fp,"%14.7e %14.7e %14.7e\n",q3[0],q3[1],q3[2]);
    fprintf(fp,"%14.7e %14.7e %14.7e\n",q4[0],q4[1],q4[2]);
    fprintf(fp,"%14.7e %14.7e %14.7e\n",q1[0],q1[1],q1[2]);
    fclose(fp);
    sprintf(buff,"facet_ti.dat");
    fp=fopen(buff,"w");
    for (i=0; i < nti; i++)
      fprintf(fp,"%14.7e %14.7e %14.7e\n",ti[i][0],ti[i][1],ti[i][2]);
    fclose(fp);
    sprintf(buff,"facet_qi.dat");
    fp=fopen(buff,"w");
    for (i=0; i < nqi; i++)
      fprintf(fp,"%14.7e %14.7e %14.7e\n",qi[i][0],qi[i][1],qi[i][2]);
    fclose(fp);
    exit(0);
  }

  // for single save nodes copy for sorting logic
  if (nti == 1)
    ti[nti++] = ti[0];
  if (nqi == 1)
    qi[nqi++] = qi[0];

  // sort the save points in each facet according to line direction
  Point swap;
  v1 = Vector(ti[0],ti[1]);
  if (v1 * line < 0.0)
  {
    swap = ti[0];
    ti[0] = ti[1];
    ti[1] = swap;
    v1 = Vector(ti[0],ti[1]);
  }
  v2 = Vector(qi[0],qi[1]);
  if (v2 * line < 0.0)
  {
    swap = qi[0];
    qi[0] = qi[1];
    qi[1] = swap;
    v2 = Vector(qi[0],qi[1]);
  }

  if ((line * Vector(lpt,ti[0])) > (line * Vector(lpt,qi[0])))
    pa = ti[0];
  else
    pa = qi[0];

  if ((line * Vector(lpt,ti[1])) < (line * Vector(lpt,qi[1])))
    pb = ti[1];
  else
    pb = qi[1];

  if ((line * Vector(lpt,pb)) < (line * Vector(lpt,pa)))
    return(0);
  else
  {
    // ensure segment is in direction of line vector
    v1 = Vector(pa,pb);
    if (v1 * line < 0.0)
    {
      swap = pa;
      pa = pb;
      pb = swap;
    }
#ifdef _DEBUG
//    sprintf(buff,"facet_segment.dat");
//    fp=fopen(buff,"w");
//    fprintf(fp,"%14.7e %14.7e %14.7e\n",pa[0],pa[1],pa[2]);
//    fprintf(fp,"%14.7e %14.7e %14.7e\n",pb[0],pb[1],pb[2]);
//    fclose(fp);
#endif
    return(1);
  }
  
}

int line_facet_intersect(Point pa, Point pb, Point p1, Point p2, Point p3, Point &p,
                         double tol)
{
  Point cg;
  Vector v1, v2, v3, norm, area;
  double da, db, rdnm, a1, a2, a3;
  double ptol = 1.0e-15;
  double tol1, tol2, tol3;

  tol1 = tol2 = tol3 = tol;
  //tol1 = distance(p2,p3)*tol*0.5;
  //tol2 = distance(p1,p3)*tol*0.5;
  //tol3 = distance(p1,p2)*tol*0.5;
  v1 = Vector(p1,p2);
  v2 = Vector(p1,p3);
  norm = (v1 % v2)*0.5;
  //double mag = MAX(1.0e-15,norm.magnitude());
  //norm /= mag;
  norm.normalize();
  cg = (p1+p2+p3)/3.0;
  v1 = Vector(cg,pa);
  v2 = Vector(cg,pb);
  da = norm * v1;
  db = norm * v2;
  v3 = Vector(pa,pb);
  v3.normalize();
  if (fabs(norm * v3) > ptol && MIN(da,db) <= 0.0 && MAX(da,db) >= 0.0)
  {
    rdnm = 1.0/MAX(1.0e-15,fabs(da)+fabs(db));
    da = fabs(da)*rdnm;
    db = 1.0-da;
    p = pa*db + pb*da; // location of possible intersection
    v1 = Vector(p,p1);
    v2 = Vector(p,p2);
    v3 = Vector(p,p3);
    area = (v2 % v3)*0.5;
    a1 = area*norm;
    area = (v3 % v1)*0.5;
    a2 = area*norm;
    area = (v1 % v2)*0.5;
    a3 = area*norm;
    if (a1 >= -tol1 && a2 >= -tol2 && a3 >= -tol3)
      return(1);
    else
      return(0);
  }

  return(0);
}

double triangle_area(Point p1, Point p2, Point p3)
{
    Vector v1 = Vector(p1,p2);
    Vector v2 = Vector(p1,p3);
    Vector norm = (v1 % v2)*0.5;
    double atotal = norm.magnitude();

  return(atotal);
}

double quadrilateral_area(Point p1, Point p2, Point p3, Point p4)
{
    Vector v1 = Vector(p1,p2);
    Vector v2 = Vector(p1,p4);
    Vector norm = (v1 % v2)*0.5;
    double atotal = norm.magnitude();
    v1 = Vector(p3,p4);
    v2 = Vector(p3,p2);
    norm = (v1 % v2)*0.5;
    atotal += norm.magnitude();

  return(atotal);
}

double triangle_aspect_ratio(Point p1, Point p2, Point p3,
                        Point &ci, Point &cc)
{
  double tol = 1.0e-15;
  double asp;

//
// compute transformation to u,v,w and back to x,y,z
//
  Vector u, v, w;  // Transformation vectors!

  u = Vector(p1,p2);
  v = Vector(p1,p3);
  w = u % v;

  u.normalize();
  w.normalize();
  v = w % u;

  Vector x, y, z;  // Transformation vectors!
  double J;
  J = (u[0]*v[1]*w[2]+u[1]*v[2]*w[0]+u[2]*v[0]*w[1]-
       u[2]*v[1]*w[0]-u[1]*v[0]*w[2]-u[0]*v[2]*w[1]);
  if (fabs(J) < tol)
  {
    ci = (p1+p2+p3)/3.0;
    cc = ci;
    asp = 1.0;
    return(asp);
  }
  J = 1.0/J;

  x[0] = J*(v[1]*w[2]-v[2]*w[1]);
  x[1] = J*(u[2]*w[1]-u[1]*w[2]);
  x[2] = J*(u[1]*v[2]-u[2]*v[1]);
  y[0] = J*(v[2]*w[0]-v[0]*w[2]);
  y[1] = J*(u[0]*w[2]-u[2]*w[0]);
  y[2] = J*(u[2]*v[0]-u[0]*v[2]);
  z[0] = J*(v[0]*w[1]-v[1]*w[0]);
  z[1] = J*(u[1]*w[0]-u[0]*w[1]);
  z[2] = J*(u[0]*v[1]-u[1]*v[0]);

  Point t1, t2, t3, tc, ti;

  t1[0] = u[0]*p1[0]+u[1]*p1[1]+u[2]*p1[2];
  t1[1] = v[0]*p1[0]+v[1]*p1[1]+v[2]*p1[2];
  t1[2] = w[0]*p1[0]+w[1]*p1[1]+w[2]*p1[2];
  t2[0] = u[0]*p2[0]+u[1]*p2[1]+u[2]*p2[2];
  t2[1] = v[0]*p2[0]+v[1]*p2[1]+v[2]*p2[2];
  t2[2] = w[0]*p2[0]+w[1]*p2[1]+w[2]*p2[2];
  t3[0] = u[0]*p3[0]+u[1]*p3[1]+u[2]*p3[2];
  t3[1] = v[0]*p3[0]+v[1]*p3[1]+v[2]*p3[2];
  t3[2] = w[0]*p3[0]+w[1]*p3[1]+w[2]*p3[2];

  double wt = (t1[2]+t2[2]+t3[2])/3.0;

//
//  compute circumscribing radius
//

  Point zero = Point(0.0,0.0,0.0);
  Vector v1 = Vector(zero,t1);
  Vector v2 = Vector(zero,t2);
  Vector v3 = Vector(zero,t3);
  double d1 = v1 * v1;
  double d2 = v2 * v2;
  double d3 = v3 * v3;
  v1 = Vector(t1,t2);
  v2 = Vector(t1,t3);
  J=v1[0]*v2[1]-v1[1]*v2[0];
  if (fabs(J) < tol)
    cc = (p1+p2+p3)/3.0;
  else
  {
    J=1.0/J;
    tc[0]=0.5*J*(v2[1]*(d2-d1)-v1[1]*(d3-d1));
    tc[1]=0.5*J*(v1[0]*(d3-d1)-v2[0]*(d2-d1));
    tc[2]=wt;
    cc[0]=tc[0]*x[0]+tc[1]*x[1]+tc[2]*x[2];
    cc[1]=tc[0]*y[0]+tc[1]*y[1]+tc[2]*y[2];
    cc[2]=tc[0]*z[0]+tc[1]*z[1]+tc[2]*z[2];
  }
  Vector v4 = Vector(cc,p1);
  double rc = v4.magnitude();

//
//compute inscribed radius
//
  Vector a, b, c;

  a[0] = t1[1]-t2[1];
  a[1] = t2[0]-t1[0];
  a[2] = 0.0;
  a.normalize();
  b[0] = t2[1]-t3[1];
  b[1] = t3[0]-t2[0];
  b[2] = 0.0;
  b.normalize();
  c[0] = t3[1]-t1[1];
  c[1] = t1[0]-t3[0];
  c[2] = 0.0;
  c.normalize();

  double ri=0.0;

  J = (-a[0]*b[1]-a[1]*c[0]-b[0]*c[1]+
        b[1]*c[0]+a[1]*b[0]+a[0]*c[1]);
  if (fabs(J) < tol)
    ci = (p1+p2+p3)/3.0;
  else
  {
    J=1.0/J;
    double a1 = a[0]*t1[0]+a[1]*t1[1];
    double b2 = b[0]*t2[0]+b[1]*t2[1];
    double c3 = c[0]*t3[0]+c[1]*t3[1];
    ri   =J*(a[0]*b[1]*c3+a[1]*b2*c[0]+a1*b[0]*c[1]-
             a1*b[1]*c[0]-a[1]*b[0]*c3-a[0]*b2*c[1]);
    ti[0]=J*(-a1*b[1]-a[1]*c3-b2*c[1]+
              b[1]*c3+a[1]*b2+a1*c[1]);
    ti[1]=J*(-a[0]*b2-a1*c[0]-b[0]*c3+
              b2*c[0]+a1*b[0]+a[0]*c3);
    ti[2]=wt;
    ci[0]=ti[0]*x[0]+ti[1]*x[1]+ti[2]*x[2];
    ci[1]=ti[0]*y[0]+ti[1]*y[1]+ti[2]*y[2];
    ci[2]=ti[0]*z[0]+ti[1]*z[1]+ti[2]*z[2];
  }

//
//compute ascpect ratio
//
  if (fabs(ri) < tol)
    asp=1.0;
  else
    asp=rc/ri/2.0;

  return(asp);

}

double quadrilateral_aspect_ratio(Point p1, Point p2, Point p3, Point p4)
{

  Point ci, cc;
  double ar;

  ar = triangle_aspect_ratio(p4,p1,p2,ci,cc);
  ar = MAX(ar,triangle_aspect_ratio(p1,p2,p3,ci,cc));
  ar = MAX(ar,triangle_aspect_ratio(p2,p3,p4,ci,cc));
  ar = MAX(ar,triangle_aspect_ratio(p3,p4,p1,ci,cc));

  return(ar);
}

double tetrahedral_volume(Point p1, Point p2, Point p3, Point p4)
{
  Vector v1 = Vector(p1,p2);
  Vector v2 = Vector(p1,p3);
  Vector v3 = Vector(p1,p4);
  double vol = scalar_triple_product(v1,v2,v3)/6.0;

  return(vol);
}

double condition_number(Point p0, Point p1, Point p2, Point p3, double w[3][3])
{
  double bn, cn;
  double a[3][3], ainv[3][3], winv[3][3], bm[3][3], cm[3][3];
  Vector e1, e2, e3;
  int flag;

  e1 = Vector(p0,p1);
  e2 = Vector(p0,p2);
  e3 = Vector(p0,p3);

  // compute condition number for tetrahedra
  ainv[0][0]=a[0][0]=e1[0];
  ainv[1][0]=a[1][0]=e1[1];
  ainv[2][0]=a[2][0]=e1[2];
  ainv[0][1]=a[0][1]=e2[0];
  ainv[1][1]=a[1][1]=e2[1];
  ainv[2][1]=a[2][1]=e2[2];
  ainv[0][2]=a[0][2]=e3[0];
  ainv[1][2]=a[1][2]=e3[1];
  ainv[2][2]=a[2][2]=e3[2];
  flag = invert3X3(ainv);

  winv[0][0]=w[0][0];
  winv[1][0]=w[1][0];
  winv[2][0]=w[2][0];
  winv[0][1]=w[0][1];
  winv[1][1]=w[1][1];
  winv[2][1]=w[2][1];
  winv[0][2]=w[0][2];
  winv[1][2]=w[1][2];
  winv[2][2]=w[2][2];
  invert3X3(winv);

  bm[0][0]=a[0][0]*winv[0][0]+a[0][1]*winv[1][0]+a[0][2]*winv[2][0];
  bm[0][1]=a[0][0]*winv[0][1]+a[0][1]*winv[1][1]+a[0][2]*winv[2][1];
  bm[0][2]=a[0][0]*winv[0][2]+a[0][1]*winv[1][2]+a[0][2]*winv[2][2];
  bm[1][0]=a[1][0]*winv[0][0]+a[1][1]*winv[1][0]+a[1][2]*winv[2][0];
  bm[1][1]=a[1][0]*winv[0][1]+a[1][1]*winv[1][1]+a[1][2]*winv[2][1];
  bm[1][2]=a[1][0]*winv[0][2]+a[1][1]*winv[1][2]+a[1][2]*winv[2][2];
  bm[2][0]=a[2][0]*winv[0][0]+a[2][1]*winv[1][0]+a[2][2]*winv[2][0];
  bm[2][1]=a[2][0]*winv[0][1]+a[2][1]*winv[1][1]+a[2][2]*winv[2][1];
  bm[2][2]=a[2][0]*winv[0][2]+a[2][1]*winv[1][2]+a[2][2]*winv[2][2];
  cm[0][0]=w[0][0]*ainv[0][0]+w[0][1]*ainv[1][0]+w[0][2]*ainv[2][0];
  cm[0][1]=w[0][0]*ainv[0][1]+w[0][1]*ainv[1][1]+w[0][2]*ainv[2][1];
  cm[0][2]=w[0][0]*ainv[0][2]+w[0][1]*ainv[1][2]+w[0][2]*ainv[2][2];
  cm[1][0]=w[1][0]*ainv[0][0]+w[1][1]*ainv[1][0]+w[1][2]*ainv[2][0];
  cm[1][1]=w[1][0]*ainv[0][1]+w[1][1]*ainv[1][1]+w[1][2]*ainv[2][1];
  cm[1][2]=w[1][0]*ainv[0][2]+w[1][1]*ainv[1][2]+w[1][2]*ainv[2][2];
  cm[2][0]=w[2][0]*ainv[0][0]+w[2][1]*ainv[1][0]+w[2][2]*ainv[2][0];
  cm[2][1]=w[2][0]*ainv[0][1]+w[2][1]*ainv[1][1]+w[2][2]*ainv[2][1];
  cm[2][2]=w[2][0]*ainv[0][2]+w[2][1]*ainv[1][2]+w[2][2]*ainv[2][2];
  bn = sqrt(bm[0][0]*bm[0][0]+bm[1][0]*bm[1][0]+bm[2][0]*bm[2][0]+
            bm[0][1]*bm[0][1]+bm[1][1]*bm[1][1]+bm[2][1]*bm[2][1]+
            bm[0][2]*bm[0][2]+bm[1][2]*bm[1][2]+bm[2][2]*bm[2][2]);
  cn = sqrt(cm[0][0]*cm[0][0]+cm[1][0]*cm[1][0]+cm[2][0]*cm[2][0]+
            cm[0][1]*cm[0][1]+cm[1][1]*cm[1][1]+cm[2][1]*cm[2][1]+
            cm[0][2]*cm[0][2]+cm[1][2]*cm[1][2]+cm[2][2]*cm[2][2]);
  cn = bn*cn/3.0;

  if (flag)
    return(cn);
  else
    return(1.0e20);
}

double tetrahedral_condition_number(Point p0, Point p1, Point p2, Point p3, int wgt)
{
  double bn, cn;
  double a[3][3], ainv[3][3], w[3][3], winv[3][3], bm[3][3], cm[3][3];
  Vector e1, e2, e3;
  int flag;

  e1 = Vector(p0,p1);
  e2 = Vector(p0,p2);
  e3 = Vector(p0,p3);

  // compute condition number for tetrahedra
  ainv[0][0]=a[0][0]=e1[0];
  ainv[1][0]=a[1][0]=e1[1];
  ainv[2][0]=a[2][0]=e1[2];
  ainv[0][1]=a[0][1]=e2[0];
  ainv[1][1]=a[1][1]=e2[1];
  ainv[2][1]=a[2][1]=e2[2];
  ainv[0][2]=a[0][2]=e3[0];
  ainv[1][2]=a[1][2]=e3[1];
  ainv[2][2]=a[2][2]=e3[2];
  flag = invert3X3(ainv);

  if (wgt)
  {
    winv[0][0]=w[0][0]=1.0;
    winv[1][0]=w[1][0]=0.0;
    winv[2][0]=w[2][0]=0.0;
    winv[0][1]=w[0][1]=0.5;
    winv[1][1]=w[1][1]=sqrt(3.0)/2.0;
    winv[2][1]=w[2][1]=0.0;
    winv[0][2]=w[0][2]=0.5;
    winv[1][2]=w[1][2]=sqrt(3.0)/6.0;
    winv[2][2]=w[2][2]=sqrt(2.0)/sqrt(3.0);
    invert3X3(winv);
    bm[0][0]=a[0][0]*winv[0][0]+a[0][1]*winv[1][0]+a[0][2]*winv[2][0];
    bm[0][1]=a[0][0]*winv[0][1]+a[0][1]*winv[1][1]+a[0][2]*winv[2][1];
    bm[0][2]=a[0][0]*winv[0][2]+a[0][1]*winv[1][2]+a[0][2]*winv[2][2];
    bm[1][0]=a[1][0]*winv[0][0]+a[1][1]*winv[1][0]+a[1][2]*winv[2][0];
    bm[1][1]=a[1][0]*winv[0][1]+a[1][1]*winv[1][1]+a[1][2]*winv[2][1];
    bm[1][2]=a[1][0]*winv[0][2]+a[1][1]*winv[1][2]+a[1][2]*winv[2][2];
    bm[2][0]=a[2][0]*winv[0][0]+a[2][1]*winv[1][0]+a[2][2]*winv[2][0];
    bm[2][1]=a[2][0]*winv[0][1]+a[2][1]*winv[1][1]+a[2][2]*winv[2][1];
    bm[2][2]=a[2][0]*winv[0][2]+a[2][1]*winv[1][2]+a[2][2]*winv[2][2];
    cm[0][0]=w[0][0]*ainv[0][0]+w[0][1]*ainv[1][0]+w[0][2]*ainv[2][0];
    cm[0][1]=w[0][0]*ainv[0][1]+w[0][1]*ainv[1][1]+w[0][2]*ainv[2][1];
    cm[0][2]=w[0][0]*ainv[0][2]+w[0][1]*ainv[1][2]+w[0][2]*ainv[2][2];
    cm[1][0]=w[1][0]*ainv[0][0]+w[1][1]*ainv[1][0]+w[1][2]*ainv[2][0];
    cm[1][1]=w[1][0]*ainv[0][1]+w[1][1]*ainv[1][1]+w[1][2]*ainv[2][1];
    cm[1][2]=w[1][0]*ainv[0][2]+w[1][1]*ainv[1][2]+w[1][2]*ainv[2][2];
    cm[2][0]=w[2][0]*ainv[0][0]+w[2][1]*ainv[1][0]+w[2][2]*ainv[2][0];
    cm[2][1]=w[2][0]*ainv[0][1]+w[2][1]*ainv[1][1]+w[2][2]*ainv[2][1];
    cm[2][2]=w[2][0]*ainv[0][2]+w[2][1]*ainv[1][2]+w[2][2]*ainv[2][2];
    bn = sqrt(bm[0][0]*bm[0][0]+bm[1][0]*bm[1][0]+bm[2][0]*bm[2][0]+
              bm[0][1]*bm[0][1]+bm[1][1]*bm[1][1]+bm[2][1]*bm[2][1]+
              bm[0][2]*bm[0][2]+bm[1][2]*bm[1][2]+bm[2][2]*bm[2][2]);
    cn = sqrt(cm[0][0]*cm[0][0]+cm[1][0]*cm[1][0]+cm[2][0]*cm[2][0]+
              cm[0][1]*cm[0][1]+cm[1][1]*cm[1][1]+cm[2][1]*cm[2][1]+
              cm[0][2]*cm[0][2]+cm[1][2]*cm[1][2]+cm[2][2]*cm[2][2]);
  } else
  {
    bn = sqrt(a[0][0]*a[0][0]+a[1][0]*a[1][0]+a[2][0]*a[2][0]+
              a[0][1]*a[0][1]+a[1][1]*a[1][1]+a[2][1]*a[2][1]+
              a[0][2]*a[0][2]+a[1][2]*a[1][2]+a[2][2]*a[2][2]);
    cn = sqrt(ainv[0][0]*ainv[0][0]+ainv[1][0]*ainv[1][0]+ainv[2][0]*ainv[2][0]+
              ainv[0][1]*ainv[0][1]+ainv[1][1]*ainv[1][1]+ainv[2][1]*ainv[2][1]+
              ainv[0][2]*ainv[0][2]+ainv[1][2]*ainv[1][2]+ainv[2][2]*ainv[2][2]);
  }
  cn = bn*cn/3.0;

  if (flag)
    return(cn);
  else
    return(1.0e20);
}

double circumcenter(Point p1, Point p2, Point p3, Point p4, Point &cc)
{
  double m[4][4];
  double a, bx, by, bz, c;

  m[0][0] = p1[0];
  m[1][0] = p2[0];
  m[2][0] = p3[0];
  m[3][0] = p4[0];
  m[0][1] = p1[1];
  m[1][1] = p2[1];
  m[2][1] = p3[1];
  m[3][1] = p4[1];
  m[0][2] = p1[2];
  m[1][2] = p2[2];
  m[2][2] = p3[2];
  m[3][2] = p4[2];
  m[0][3] = 1.0;
  m[1][3] = 1.0;
  m[2][3] = 1.0;
  m[3][3] = 1.0;
  a = determinant_4X4(m);

  m[0][0] = p1[0]*p1[0]+p1[1]*p1[1]+p1[2]*p1[2];
  m[1][0] = p2[0]*p2[0]+p2[1]*p2[1]+p2[2]*p2[2];
  m[2][0] = p3[0]*p3[0]+p3[1]*p3[1]+p3[2]*p3[2];
  m[3][0] = p4[0]*p4[0]+p4[1]*p4[1]+p4[2]*p4[2];
  bx = determinant_4X4(m);

  m[0][1] = p1[0];
  m[1][1] = p2[0];
  m[2][1] = p3[0];
  m[3][1] = p4[0];
  by = determinant_4X4(m);

  m[0][2] = p1[1];
  m[1][2] = p2[1];
  m[2][2] = p3[1];
  m[3][2] = p4[1];
  bz = determinant_4X4(m);

  m[0][3] = p1[2];
  m[1][3] = p2[2];
  m[2][3] = p3[2];
  m[3][3] = p4[2];
  c = determinant_4X4(m);

  double r;
  if (fabs(a) < 1.0e-20)
  {
    cc = (p1+p2+p3+p4)*0.25;
    r = 1.0e20;
    return(r);
  }

  cc = Point(0.5*bx/a,-0.5*by/a,0.5*bz/a);

  r = (bx*bx+by*by+bz*bz-4.0*a*c);
  if (r < 0.0)
    r = 1.0e20;
  else
    r = 0.5*sqrt(r)/fabs(a);

  return(r);
}

double tetrahedral_aspect_ratio(Point p1, Point p2, Point p3, Point p4,
                        Point &ci, Point &cc)
{
  double tol = 1.0e-20;
//
//  compute circumscribing radius
//
  double m[4][4];
  double a, bx, by, bz, c;
  double mag, rc;

  m[0][0] = p1[0];
  m[1][0] = p2[0];
  m[2][0] = p3[0];
  m[3][0] = p4[0];
  m[0][1] = p1[1];
  m[1][1] = p2[1];
  m[2][1] = p3[1];
  m[3][1] = p4[1];
  m[0][2] = p1[2];
  m[1][2] = p2[2];
  m[2][2] = p3[2];
  m[3][2] = p4[2];
  m[0][3] = 1.0;
  m[1][3] = 1.0;
  m[2][3] = 1.0;
  m[3][3] = 1.0;
  a = determinant_4X4(m);

  m[0][0] = p1[0]*p1[0]+p1[1]*p1[1]+p1[2]*p1[2];
  m[1][0] = p2[0]*p2[0]+p2[1]*p2[1]+p2[2]*p2[2];
  m[2][0] = p3[0]*p3[0]+p3[1]*p3[1]+p3[2]*p3[2];
  m[3][0] = p4[0]*p4[0]+p4[1]*p4[1]+p4[2]*p4[2];
  bx = determinant_4X4(m);

  m[0][1] = p1[0];
  m[1][1] = p2[0];
  m[2][1] = p3[0];
  m[3][1] = p4[0];
  by = determinant_4X4(m);

  m[0][2] = p1[1];
  m[1][2] = p2[1];
  m[2][2] = p3[1];
  m[3][2] = p4[1];
  bz = determinant_4X4(m);

  m[0][3] = p1[2];
  m[1][3] = p2[2];
  m[2][3] = p3[2];
  m[3][3] = p4[2];
  c = determinant_4X4(m);

  if (fabs(a) > tol)
    cc = Point(0.5*bx/a,-0.5*by/a,0.5*bz/a);
  else
    cc = (p1+p2+p3+p4)*0.25;

  rc = (distance(cc,p1)+distance(cc,p2)+
        distance(cc,p3)+distance(cc,p4))*0.25;
  //rc = 0.5*sqrt(bx*bx+by*by+bz*bz-4.0*a*c)/fabs(a);

//
//compute inscribed radius
//
  Vector v1 = Vector(p1,p2);
  Vector v2 = Vector(p1,p3);
  Vector v3 = Vector(p1,p4);
  Vector v4 = Vector(p2,p3);
  Vector v5 = Vector(p2,p4);
  double a1=(v5[1]*v4[2]-v5[2]*v4[1]);
  double b1=(v5[2]*v4[0]-v5[0]*v4[2]);
  double c1=(v5[0]*v4[1]-v5[1]*v4[0]);
  mag = 1.0/MAX(tol,sqrt(a1*a1+b1*b1+c1*c1));
  a1 *= mag;
  b1 *= mag;
  c1 *= mag;
  double a2=(v2[1]*v3[2]-v2[2]*v3[1]);
  double b2=(v2[2]*v3[0]-v2[0]*v3[2]);
  double c2=(v2[0]*v3[1]-v2[1]*v3[0]);
  mag = 1.0/MAX(tol,sqrt(a2*a2+b2*b2+c2*c2));
  a2 *= mag;
  b2 *= mag;
  c2 *= mag;
  double a3=(v3[1]*v1[2]-v3[2]*v1[1]);
  double b3=(v3[2]*v1[0]-v3[0]*v1[2]);
  double c3=(v3[0]*v1[1]-v3[1]*v1[0]);
  mag = 1.0/MAX(tol,sqrt(a3*a3+b3*b3+c3*c3));
  a3 *= mag;
  b3 *= mag;
  c3 *= mag;
  double a4=(v1[1]*v2[2]-v1[2]*v2[1]);
  double b4=(v1[2]*v2[0]-v1[0]*v2[2]);
  double c4=(v1[0]*v2[1]-v1[1]*v2[0]);
  mag = 1.0/MAX(tol,sqrt(a4*a4+b4*b4+c4*c4));
  a4 *= mag;
  b4 *= mag;
  c4 *= mag;
  double d1= p4[0]*a1+p4[1]*b1+p4[2]*c1;
  double d2= p3[0]*a2+p3[1]*b2+p3[2]*c2;
  double d3= p2[0]*a3+p2[1]*b3+p2[2]*c3;
  double d4= p1[0]*a4+p1[1]*b4+p1[2]*c4;
  double et=a2*b3*c4+b2*c3*a4+c2*a3*b4-c2*b3*a4-b2*a3*c4-a2*c3*b4;
  double ft=a1*b3*c4+b1*c3*a4+c1*a3*b4-c1*b3*a4-b1*a3*c4-a1*c3*b4;
  double gt=a1*b2*c4+b1*c2*a4+c1*a2*b4-c1*b2*a4-b1*a2*c4-a1*c2*b4;
  double ht=a1*b2*c3+b1*c2*a3+c1*a2*b3-c1*b2*a3-b1*a2*c3-a1*c2*b3;
  double rd=(et-ft+gt-ht);
  double ri;
  if (fabs(rd) < tol)
    ri = tol;
  else
    ri=((-d1*et+d2*ft-d3*gt+d4*ht)/rd);
//
//compute ascpect ratio
//
  double tet_aspect;
  if (fabs(ri) < tol)
    tet_aspect=1.0/tol;
  else
    tet_aspect=rc/ri/3.0;

  if (fabs(rd) < tol || fabs(ri) < tol)
    ci = (p1+p2+p3+p4)*0.25;
  else
  {
    rd=1.0/rd;
    et=d2*b3*c4+b2*c3*d4+c2*d3*b4-c2*b3*d4-b2*d3*c4-d2*c3*b4;
    ft=d1*b3*c4+b1*c3*d4+c1*d3*b4-c1*b3*d4-b1*d3*c4-d1*c3*b4;
    gt=d1*b2*c4+b1*c2*d4+c1*d2*b4-c1*b2*d4-b1*d2*c4-d1*c2*b4;
    ht=d1*b2*c3+b1*c2*d3+c1*d2*b3-c1*b2*d3-b1*d2*c3-d1*c2*b3;
    ci[0]=(et-ft+gt-ht)*rd;
    et=a2*d3*c4+d2*c3*a4+c2*a3*d4-c2*d3*a4-d2*a3*c4-a2*c3*d4;
    ft=a1*d3*c4+d1*c3*a4+c1*a3*d4-c1*d3*a4-d1*a3*c4-a1*c3*d4;
    gt=a1*d2*c4+d1*c2*a4+c1*a2*d4-c1*d2*a4-d1*a2*c4-a1*c2*d4;
    ht=a1*d2*c3+d1*c2*a3+c1*a2*d3-c1*d2*a3-d1*a2*c3-a1*c2*d3;
    ci[1]=(et-ft+gt-ht)*rd;
    et=a2*b3*d4+b2*d3*a4+d2*a3*b4-d2*b3*a4-b2*a3*d4-a2*d3*b4;
    ft=a1*b3*d4+b1*d3*a4+d1*a3*b4-d1*b3*a4-b1*a3*d4-a1*d3*b4;
    gt=a1*b2*d4+b1*d2*a4+d1*a2*b4-d1*b2*a4-b1*a2*d4-a1*d2*b4;
    ht=a1*b2*d3+b1*d2*a3+d1*a2*b3-d1*b2*a3-b1*a2*d3-a1*d2*b3;
    ci[2]=(et-ft+gt-ht)*rd;
  }

  return(tet_aspect);

}

double tetrahedral_aspect_ratio(Point p1, Point p2, Point p3, Point p4)
{
  double a1, a2, a3, a4, a, b, c;
  double vol, e1, e2, e3, e4, e5, e6;
  Vector v1, v2, v3;

  vol = tetrahedral_volume(p1, p2, p3, p4);
  e1 = distance(p1,p2);
  e2 = distance(p1,p3);
  e3 = distance(p1,p4);
  e4 = distance(p2,p3);
  e5 = distance(p2,p4);
  e6 = distance(p3,p4);
  v1 = Vector(p1,p2);
  v2 = Vector(p1,p3);
  v3 = v1 % v2;
  a1 = v3.magnitude()*0.5;
  v1 = Vector(p1,p2);
  v2 = Vector(p1,p4);
  v3 = v1 % v2;
  a2 = v3.magnitude()*0.5;
  v1 = Vector(p1,p3);
  v2 = Vector(p1,p4);
  v3 = v1 % v2;
  a3 = v3.magnitude()*0.5;
  v1 = Vector(p2,p3);
  v2 = Vector(p2,p4);
  v3 = v1 % v2;
  a4 = v3.magnitude()*0.5;

  a = e1*e6;
  b = e2*e5;
  c = e3*e4;
  double rin = 3.0*vol/MAX(1.0e-15,a1+a2+a3+a4);
  double rout= sqrt((a+b+c)*(a+b-c)*(a+c-b)*(b+c-a))/24.0/MAX(1.0e-15,vol);
  double ar = rout/MAX(1.0e-15,rin)/3.0;

  return(ar);
}

double pyramid_aspect_ratio(Point p1, Point p2, Point p3, Point p4, Point p5)
{
  Point fc = (p1 + p2 + p3 + p4)*0.25;
  double L, W, H, ar, mn, mx;
  H = distance(p5,fc);
  L = distance((p1+p4)*0.5,(p2+p3)*0.5);
  W = distance((p1+p2)*0.5,(p3+p4)*0.5);
  mx = MAX(H,MAX(L,W));
  mn = MIN(H,MIN(L,W));
  ar = mx/MAX(1.0e-20,mn);
  return(ar);
}

double prism_aspect_ratio(Point p1, Point p2, Point p3, Point p4, Point p5, Point p6)
{
  Point ci, cc;
  double ar, H, R, mn, mx;

  ar = triangle_aspect_ratio((p1+p4)*0.5,(p2+p5)*0.5,(p3+p6)*0.5,ci,cc);
  R = distance(cc,(p1+p4)*0.5);
  H = distance((p1+p2+p3)/3.0,(p4+p5+p6)/3.0);
  mn = MIN(R,H);
  mx = MAX(R,H);
  ar = mx/MAX(1.0e-20,mn);
  return(ar);
}

double hexahedral_aspect_ratio(Point p1, Point p2, Point p3, Point p4, Point p5, Point p6, Point p7, Point p8)
{
  double L, W, H, ar, mn, mx;

  W = distance((p1+p4+p5+p8)*0.25,(p2+p3+p6+p7)*0.25);
  L = distance((p1+p2+p5+p6)*0.25,(p3+p4+p7+p8)*0.25);
  H = distance((p1+p2+p3+p4)*0.25,(p5+p6+p7+p8)*0.25);
  mx = MAX(H,MAX(L,W));
  mn = MIN(H,MIN(L,W));
  ar = mx/MAX(1.0e-20,mn);
  return(ar);
}

double pyramid_volume(Point p1, Point p2, Point p3, Point p4, Point p5)
{
  Point fc = (p1 + p2 + p3 + p4)*0.25;
  Vector v1, v2, v3;
  double vol;
  v1 = Vector(fc,p1);
  v2 = Vector(fc,p2);
  v3 = Vector(fc,p5);
  vol = scalar_triple_product(v1,v2,v3)/6.0;
  v1 = Vector(fc,p2);
  v2 = Vector(fc,p3);
  v3 = Vector(fc,p5);
  vol += scalar_triple_product(v1,v2,v3)/6.0;
  v1 = Vector(fc,p3);
  v2 = Vector(fc,p4);
  v3 = Vector(fc,p5);
  vol += scalar_triple_product(v1,v2,v3)/6.0;
  v1 = Vector(fc,p4);
  v2 = Vector(fc,p1);
  v3 = Vector(fc,p5);
  vol += scalar_triple_product(v1,v2,v3)/6.0;

  return(vol);
}

double prism_volume(Point p1, Point p2, Point p3, Point p4, Point p5, Point p6)
{
  Point cg = (p1 + p2 + p3 + p4 + p5 + p6)/6.0;
  Point fc;
  Vector v1, v2, v3;
  double vol = 0.0;
  // side 1
  fc = (p1 + p2 + p4 + p5)*0.25;
  v1 = Vector(fc,p2);
  v2 = Vector(fc,p1);
  v3 = Vector(fc,cg);
  vol += scalar_triple_product(v1,v2,v3)/6.0;
  v1 = Vector(fc,p1);
  v2 = Vector(fc,p4);
  v3 = Vector(fc,cg);
  vol += scalar_triple_product(v1,v2,v3)/6.0;
  v1 = Vector(fc,p4);
  v2 = Vector(fc,p5);
  v3 = Vector(fc,cg);
  vol += scalar_triple_product(v1,v2,v3)/6.0;
  v1 = Vector(fc,p5);
  v2 = Vector(fc,p2);
  v3 = Vector(fc,cg);
  vol += scalar_triple_product(v1,v2,v3)/6.0;
  // side 2
  fc = (p2 + p3 + p5 + p6)*0.25;
  v1 = Vector(fc,p3);
  v2 = Vector(fc,p2);
  v3 = Vector(fc,cg);
  vol += scalar_triple_product(v1,v2,v3)/6.0;
  v1 = Vector(fc,p2);
  v2 = Vector(fc,p5);
  v3 = Vector(fc,cg);
  vol += scalar_triple_product(v1,v2,v3)/6.0;
  v1 = Vector(fc,p5);
  v2 = Vector(fc,p6);
  v3 = Vector(fc,cg);
  vol += scalar_triple_product(v1,v2,v3)/6.0;
  v1 = Vector(fc,p6);
  v2 = Vector(fc,p3);
  v3 = Vector(fc,cg);
  vol += scalar_triple_product(v1,v2,v3)/6.0;
  // side 3
  fc = (p3 + p1 + p6 + p4)*0.25;
  v1 = Vector(fc,p1);
  v2 = Vector(fc,p3);
  v3 = Vector(fc,cg);
  vol += scalar_triple_product(v1,v2,v3)/6.0;
  v1 = Vector(fc,p3);
  v2 = Vector(fc,p6);
  v3 = Vector(fc,cg);
  vol += scalar_triple_product(v1,v2,v3)/6.0;
  v1 = Vector(fc,p6);
  v2 = Vector(fc,p4);
  v3 = Vector(fc,cg);
  vol += scalar_triple_product(v1,v2,v3)/6.0;
  v1 = Vector(fc,p4);
  v2 = Vector(fc,p1);
  v3 = Vector(fc,cg);
  vol += scalar_triple_product(v1,v2,v3)/6.0;
  // side 4
  v1 = Vector(p1,p2);
  v2 = Vector(p1,p3);
  v3 = Vector(p1,cg);
  vol += scalar_triple_product(v1,v2,v3)/6.0;
  // side 5
  v1 = Vector(p4,p6);
  v2 = Vector(p4,p5);
  v3 = Vector(p4,cg);
  vol += scalar_triple_product(v1,v2,v3)/6.0;

  return(vol);
}

double hexahedral_volume(Point p1, Point p2, Point p3, Point p4,
                         Point p5, Point p6, Point p7, Point p8)
{
  Point cg = (p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8)/8.0;
  Point fc;
  Vector v1, v2, v3;
  double vol = 0.0;
  fc = (p1 + p2 + p3 + p4)*0.25;
  v1 = Vector(fc,p1);
  v2 = Vector(fc,p2);
  v3 = Vector(fc,cg);
  vol += scalar_triple_product(v1,v2,v3)/6.0;
  v1 = Vector(fc,p2);
  v2 = Vector(fc,p3);
  v3 = Vector(fc,cg);
  vol += scalar_triple_product(v1,v2,v3)/6.0;
  v1 = Vector(fc,p3);
  v2 = Vector(fc,p4);
  v3 = Vector(fc,cg);
  vol += scalar_triple_product(v1,v2,v3)/6.0;
  v1 = Vector(fc,p4);
  v2 = Vector(fc,p1);
  v3 = Vector(fc,cg);
  vol += scalar_triple_product(v1,v2,v3)/6.0;
  fc = (p1 + p2 + p5 + p6)*0.25;
  v1 = Vector(fc,p1);
  v2 = Vector(fc,p5);
  v3 = Vector(fc,cg);
  vol += scalar_triple_product(v1,v2,v3)/6.0;
  v1 = Vector(fc,p5);
  v2 = Vector(fc,p6);
  v3 = Vector(fc,cg);
  vol += scalar_triple_product(v1,v2,v3)/6.0;
  v1 = Vector(fc,p6);
  v2 = Vector(fc,p2);
  v3 = Vector(fc,cg);
  vol += scalar_triple_product(v1,v2,v3)/6.0;
  v1 = Vector(fc,p2);
  v2 = Vector(fc,p1);
  v3 = Vector(fc,cg);
  vol += scalar_triple_product(v1,v2,v3)/6.0;
  fc = (p2 + p3 + p6 + p7)*0.25;
  v1 = Vector(fc,p2);
  v2 = Vector(fc,p6);
  v3 = Vector(fc,cg);
  vol += scalar_triple_product(v1,v2,v3)/6.0;
  v1 = Vector(fc,p6);
  v2 = Vector(fc,p7);
  v3 = Vector(fc,cg);
  vol += scalar_triple_product(v1,v2,v3)/6.0;
  v1 = Vector(fc,p7);
  v2 = Vector(fc,p3);
  v3 = Vector(fc,cg);
  vol += scalar_triple_product(v1,v2,v3)/6.0;
  v1 = Vector(fc,p3);
  v2 = Vector(fc,p2);
  v3 = Vector(fc,cg);
  vol += scalar_triple_product(v1,v2,v3)/6.0;
  fc = (p3 + p4 + p7 + p8)*0.25;
  v1 = Vector(fc,p3);
  v2 = Vector(fc,p7);
  v3 = Vector(fc,cg);
  vol += scalar_triple_product(v1,v2,v3)/6.0;
  v1 = Vector(fc,p7);
  v2 = Vector(fc,p8);
  v3 = Vector(fc,cg);
  vol += scalar_triple_product(v1,v2,v3)/6.0;
  v1 = Vector(fc,p8);
  v2 = Vector(fc,p4);
  v3 = Vector(fc,cg);
  vol += scalar_triple_product(v1,v2,v3)/6.0;
  v1 = Vector(fc,p4);
  v2 = Vector(fc,p3);
  v3 = Vector(fc,cg);
  vol += scalar_triple_product(v1,v2,v3)/6.0;
  fc = (p1 + p4 + p5 + p8)*0.25;
  v1 = Vector(fc,p4);
  v2 = Vector(fc,p8);
  v3 = Vector(fc,cg);
  vol += scalar_triple_product(v1,v2,v3)/6.0;
  v1 = Vector(fc,p8);
  v2 = Vector(fc,p5);
  v3 = Vector(fc,cg);
  vol += scalar_triple_product(v1,v2,v3)/6.0;
  v1 = Vector(fc,p5);
  v2 = Vector(fc,p1);
  v3 = Vector(fc,cg);
  vol += scalar_triple_product(v1,v2,v3)/6.0;
  v1 = Vector(fc,p1);
  v2 = Vector(fc,p4);
  v3 = Vector(fc,cg);
  vol += scalar_triple_product(v1,v2,v3)/6.0;
  fc = (p5 + p6 + p7 + p8)*0.25;
  v1 = Vector(fc,p5);
  v2 = Vector(fc,p8);
  v3 = Vector(fc,cg);
  vol += scalar_triple_product(v1,v2,v3)/6.0;
  v1 = Vector(fc,p8);
  v2 = Vector(fc,p7);
  v3 = Vector(fc,cg);
  vol += scalar_triple_product(v1,v2,v3)/6.0;
  v1 = Vector(fc,p7);
  v2 = Vector(fc,p6);
  v3 = Vector(fc,cg);
  vol += scalar_triple_product(v1,v2,v3)/6.0;
  v1 = Vector(fc,p6);
  v2 = Vector(fc,p5);
  v3 = Vector(fc,cg);
  vol += scalar_triple_product(v1,v2,v3)/6.0;

  return(vol);
}

Vector tetrahedral_gradient(Point p0, Point p1, Point p2, Point p3,
                            double f0, double f1, double f2, double f3)
{
  double cx, cy, cz, vol, ff;
  Vector v1, v2, area, grad;

  vol = tetrahedral_volume(p0,p1,p2,p3);

  cx = cy = cz = 0.0;

  // face 0
  v1 = Vector(p0,p2);
  v2 = Vector(p0,p1);
  area = (v1 % v2)*0.5;
  ff = (f0+f1+f2)/3.0;
  cx += ff*area[0];
  cy += ff*area[1];
  cz += ff*area[2];
  // face 1
  v1 = Vector(p0,p1);
  v2 = Vector(p0,p3);
  area = (v1 % v2)*0.5;
  ff = (f0+f1+f3)/3.0;
  cx += ff*area[0];
  cy += ff*area[1];
  cz += ff*area[2];
  // face 2
  v1 = Vector(p1,p2);
  v2 = Vector(p1,p3);
  area = (v1 % v2)*0.5;
  ff = (f1+f2+f3)/3.0;
  cx += ff*area[0];
  cy += ff*area[1];
  cz += ff*area[2];
  // face 3
  v1 = Vector(p2,p0);
  v2 = Vector(p2,p3);
  area = (v1 % v2)*0.5;
  ff = (f2+f0+f3)/3.0;
  cx += ff*area[0];
  cy += ff*area[1];
  cz += ff*area[2];

  grad = Vector(cx/vol,cy/vol,cz/vol);

  return(grad);
}

Vector pyramid_gradient(Point p0, Point p1, Point p2, Point p3, Point p4,
                        double f0, double f1, double f2, double f3, double f4)
{
  double cx, cy, cz, vol, ff, fcg;
  Vector v1, v2, area, grad;
  Point cg;

  vol = pyramid_volume(p0,p1,p2,p3,p4);

  cx = cy = cz = 0.0;

  // face 0
  cg = (p0+p1+p2+p3)*0.25;
  fcg = (f0+f1+f2+f3)*0.25;
  v1 = Vector(cg,p0);
  v2 = Vector(cg,p3);
  area = (v1 % v2)*0.5;
  ff = (fcg+f0+f3)/3.0;
  cx += ff*area[0];
  cy += ff*area[1];
  cz += ff*area[2];
  v1 = Vector(cg,p3);
  v2 = Vector(cg,p2);
  area = (v1 % v2)*0.5;
  ff = (fcg+f3+f2)/3.0;
  cx += ff*area[0];
  cy += ff*area[1];
  cz += ff*area[2];
  v1 = Vector(cg,p2);
  v2 = Vector(cg,p1);
  area = (v1 % v2)*0.5;
  ff = (fcg+f2+f1)/3.0;
  cx += ff*area[0];
  cy += ff*area[1];
  cz += ff*area[2];
  v1 = Vector(cg,p1);
  v2 = Vector(cg,p0);
  area = (v1 % v2)*0.5;
  ff = (fcg+f1+f0)/3.0;
  cx += ff*area[0];
  cy += ff*area[1];
  cz += ff*area[2];
  // face 1
  v1 = Vector(p4,p0);
  v2 = Vector(p4,p1);
  area = (v1 % v2)*0.5;
  ff = (f0+f1+f4)/3.0;
  cx += ff*area[0];
  cy += ff*area[1];
  cz += ff*area[2];
  // face 2
  v1 = Vector(p4,p1);
  v2 = Vector(p4,p2);
  area = (v1 % v2)*0.5;
  ff = (f1+f2+f4)/3.0;
  cx += ff*area[0];
  cy += ff*area[1];
  cz += ff*area[2];
  // face 3
  v1 = Vector(p4,p2);
  v2 = Vector(p4,p3);
  area = (v1 % v2)*0.5;
  ff = (f2+f3+f4)/3.0;
  cx += ff*area[0];
  cy += ff*area[1];
  cz += ff*area[2];
  // face 4
  v1 = Vector(p4,p3);
  v2 = Vector(p4,p0);
  area = (v1 % v2)*0.5;
  ff = (f0+f3+f4)/3.0;
  cx += ff*area[0];
  cy += ff*area[1];
  cz += ff*area[2];

  grad = Vector(cx/vol,cy/vol,cz/vol);

  return(grad);
}

Vector prism_gradient(Point p0, Point p1, Point p2,
                      Point p3, Point p4, Point p5,
                      double f0, double f1, double f2,
                      double f3, double f4, double f5)
{
  double cx, cy, cz, vol, ff, fcg;
  Vector v1, v2, area, grad;
  Point cg;

  vol = prism_volume(p0,p1,p2,p3,p4,p5);

  cx = cy = cz = 0.0;

  // face 0
  cg = (p0+p1+p4+p3)*0.25;
  fcg = (f0+f1+f4+f3)*0.25;
  v1 = Vector(cg,p0);
  v2 = Vector(cg,p1);
  area = (v1 % v2)*0.5;
  ff = (fcg+f0+f1)/3.0;
  cx += ff*area[0];
  cy += ff*area[1];
  cz += ff*area[2];
  v1 = Vector(cg,p1);
  v2 = Vector(cg,p4);
  area = (v1 % v2)*0.5;
  ff = (fcg+f1+f4)/3.0;
  cx += ff*area[0];
  cy += ff*area[1];
  cz += ff*area[2];
  v1 = Vector(cg,p4);
  v2 = Vector(cg,p3);
  area = (v1 % v2)*0.5;
  ff = (fcg+f4+f3)/3.0;
  cx += ff*area[0];
  cy += ff*area[1];
  cz += ff*area[2];
  v1 = Vector(cg,p3);
  v2 = Vector(cg,p0);
  area = (v1 % v2)*0.5;
  ff = (fcg+f3+f0)/3.0;
  cx += ff*area[0];
  cy += ff*area[1];
  cz += ff*area[2];
  // face 1
  cg = (p1+p2+p5+p4)*0.25;
  fcg = (f1+f2+f5+f4)*0.25;
  v1 = Vector(cg,p1);
  v2 = Vector(cg,p2);
  area = (v1 % v2)*0.5;
  ff = (fcg+f1+f2)/3.0;
  cx += ff*area[0];
  cy += ff*area[1];
  cz += ff*area[2];
  v1 = Vector(cg,p2);
  v2 = Vector(cg,p5);
  area = (v1 % v2)*0.5;
  ff = (fcg+f2+f5)/3.0;
  cx += ff*area[0];
  cy += ff*area[1];
  cz += ff*area[2];
  v1 = Vector(cg,p5);
  v2 = Vector(cg,p4);
  area = (v1 % v2)*0.5;
  ff = (fcg+f5+f4)/3.0;
  cx += ff*area[0];
  cy += ff*area[1];
  cz += ff*area[2];
  v1 = Vector(cg,p4);
  v2 = Vector(cg,p1);
  area = (v1 % v2)*0.5;
  ff = (fcg+f4+f1)/3.0;
  cx += ff*area[0];
  cy += ff*area[1];
  cz += ff*area[2];
  // face 2
  cg = (p2+p0+p3+p5)*0.25;
  fcg = (f2+f0+f3+f5)*0.25;
  v1 = Vector(cg,p2);
  v2 = Vector(cg,p0);
  area = (v1 % v2)*0.5;
  ff = (fcg+f2+f0)/3.0;
  cx += ff*area[0];
  cy += ff*area[1];
  cz += ff*area[2];
  v1 = Vector(cg,p0);
  v2 = Vector(cg,p3);
  area = (v1 % v2)*0.5;
  ff = (fcg+f0+f3)/3.0;
  cx += ff*area[0];
  cy += ff*area[1];
  cz += ff*area[2];
  v1 = Vector(cg,p3);
  v2 = Vector(cg,p5);
  area = (v1 % v2)*0.5;
  ff = (fcg+f3+f5)/3.0;
  cx += ff*area[0];
  cy += ff*area[1];
  cz += ff*area[2];
  v1 = Vector(cg,p5);
  v2 = Vector(cg,p2);
  area = (v1 % v2)*0.5;
  ff = (fcg+f5+f2)/3.0;
  cx += ff*area[0];
  cy += ff*area[1];
  cz += ff*area[2];
  // face 3
  v1 = Vector(p0,p2);
  v2 = Vector(p0,p1);
  area = (v1 % v2)*0.5;
  ff = (f0+f2+f1)/3.0;
  cx += ff*area[0];
  cy += ff*area[1];
  cz += ff*area[2];
  // face 4
  v1 = Vector(p3,p4);
  v2 = Vector(p3,p5);
  area = (v1 % v2)*0.5;
  ff = (f3+f4+f5)/3.0;
  cx += ff*area[0];
  cy += ff*area[1];
  cz += ff*area[2];

  grad = Vector(cx/vol,cy/vol,cz/vol);

  return(grad);
}

Vector hexahedral_gradient(Point p0, Point p1, Point p2, Point p3,
                           Point p4, Point p5, Point p6, Point p7,
                           double f0, double f1, double f2, double f3,
                           double f4, double f5, double f6, double f7)
{
  double cx, cy, cz, vol, ff, fcg;
  Vector v1, v2, area, grad;
  Point cg;

  vol = hexahedral_volume(p0,p1,p2,p3,p4,p5,p6,p7);

  cx = cy = cz = 0.0;

  // face 0
  cg = (p0+p3+p2+p1)*0.25;
  fcg = (f0+f3+f2+f1)*0.25;
  v1 = Vector(cg,p0);
  v2 = Vector(cg,p3);
  area = (v1 % v2)*0.5;
  ff = (fcg+f0+f3)/3.0;
  cx += ff*area[0];
  cy += ff*area[1];
  cz += ff*area[2];
  v1 = Vector(cg,p3);
  v2 = Vector(cg,p2);
  area = (v1 % v2)*0.5;
  ff = (fcg+f3+f2)/3.0;
  cx += ff*area[0];
  cy += ff*area[1];
  cz += ff*area[2];
  v1 = Vector(cg,p2);
  v2 = Vector(cg,p1);
  area = (v1 % v2)*0.5;
  ff = (fcg+f2+f1)/3.0;
  cx += ff*area[0];
  cy += ff*area[1];
  cz += ff*area[2];
  v1 = Vector(cg,p1);
  v2 = Vector(cg,p0);
  area = (v1 % v2)*0.5;
  ff = (fcg+f1+f0)/3.0;
  cx += ff*area[0];
  cy += ff*area[1];
  cz += ff*area[2];
  // face 1
  cg = (p0+p1+p5+p4)*0.25;
  fcg = (f0+f1+f5+f4)*0.25;
  v1 = Vector(cg,p0);
  v2 = Vector(cg,p1);
  area = (v1 % v2)*0.5;
  ff = (fcg+f0+f1)/3.0;
  cx += ff*area[0];
  cy += ff*area[1];
  cz += ff*area[2];
  v1 = Vector(cg,p1);
  v2 = Vector(cg,p5);
  area = (v1 % v2)*0.5;
  ff = (fcg+f1+f5)/3.0;
  cx += ff*area[0];
  cy += ff*area[1];
  cz += ff*area[2];
  v1 = Vector(cg,p5);
  v2 = Vector(cg,p4);
  area = (v1 % v2)*0.5;
  ff = (fcg+f5+f4)/3.0;
  cx += ff*area[0];
  cy += ff*area[1];
  cz += ff*area[2];
  v1 = Vector(cg,p4);
  v2 = Vector(cg,p0);
  area = (v1 % v2)*0.5;
  ff = (fcg+f4+f0)/3.0;
  cx += ff*area[0];
  cy += ff*area[1];
  cz += ff*area[2];
  // face 2
  cg = (p1+p2+p6+p5)*0.25;
  fcg = (f1+f2+f6+f5)*0.25;
  v1 = Vector(cg,p1);
  v2 = Vector(cg,p2);
  area = (v1 % v2)*0.5;
  ff = (fcg+f1+f2)/3.0;
  cx += ff*area[0];
  cy += ff*area[1];
  cz += ff*area[2];
  v1 = Vector(cg,p2);
  v2 = Vector(cg,p6);
  area = (v1 % v2)*0.5;
  ff = (fcg+f2+f6)/3.0;
  cx += ff*area[0];
  cy += ff*area[1];
  cz += ff*area[2];
  v1 = Vector(cg,p6);
  v2 = Vector(cg,p5);
  area = (v1 % v2)*0.5;
  ff = (fcg+f6+f5)/3.0;
  cx += ff*area[0];
  cy += ff*area[1];
  cz += ff*area[2];
  v1 = Vector(cg,p5);
  v2 = Vector(cg,p1);
  area = (v1 % v2)*0.5;
  ff = (fcg+f5+f1)/3.0;
  cx += ff*area[0];
  cy += ff*area[1];
  cz += ff*area[2];
  // face 3
  cg = (p2+p3+p7+p6)*0.25;
  fcg = (f2+f3+f7+f6)*0.25;
  v1 = Vector(cg,p2);
  v2 = Vector(cg,p3);
  area = (v1 % v2)*0.5;
  ff = (fcg+f2+f3)/3.0;
  cx += ff*area[0];
  cy += ff*area[1];
  cz += ff*area[2];
  v1 = Vector(cg,p3);
  v2 = Vector(cg,p7);
  area = (v1 % v2)*0.5;
  ff = (fcg+f3+f7)/3.0;
  cx += ff*area[0];
  cy += ff*area[1];
  cz += ff*area[2];
  v1 = Vector(cg,p7);
  v2 = Vector(cg,p6);
  area = (v1 % v2)*0.5;
  ff = (fcg+f7+f6)/3.0;
  cx += ff*area[0];
  cy += ff*area[1];
  cz += ff*area[2];
  v1 = Vector(cg,p6);
  v2 = Vector(cg,p2);
  area = (v1 % v2)*0.5;
  ff = (fcg+f6+f2)/3.0;
  cx += ff*area[0];
  cy += ff*area[1];
  cz += ff*area[2];
  // face 4
  cg = (p0+p4+p7+p3)*0.25;
  fcg = (f0+f4+f7+f3)*0.25;
  v1 = Vector(cg,p0);
  v2 = Vector(cg,p4);
  area = (v1 % v2)*0.5;
  ff = (fcg+f0+f4)/3.0;
  cx += ff*area[0];
  cy += ff*area[1];
  cz += ff*area[2];
  v1 = Vector(cg,p4);
  v2 = Vector(cg,p7);
  area = (v1 % v2)*0.5;
  ff = (fcg+f4+f7)/3.0;
  cx += ff*area[0];
  cy += ff*area[1];
  cz += ff*area[2];
  v1 = Vector(cg,p7);
  v2 = Vector(cg,p3);
  area = (v1 % v2)*0.5;
  ff = (fcg+f7+f3)/3.0;
  cx += ff*area[0];
  cy += ff*area[1];
  cz += ff*area[2];
  v1 = Vector(cg,p3);
  v2 = Vector(cg,p0);
  area = (v1 % v2)*0.5;
  ff = (fcg+f3+f0)/3.0;
  cx += ff*area[0];
  cy += ff*area[1];
  cz += ff*area[2];
  // face 5
  cg = (p4+p5+p6+p7)*0.25;
  fcg = (f4+f5+f6+f7)*0.25;
  v1 = Vector(cg,p4);
  v2 = Vector(cg,p5);
  area = (v1 % v2)*0.5;
  ff = (fcg+f4+f5)/3.0;
  cx += ff*area[0];
  cy += ff*area[1];
  cz += ff*area[2];
  v1 = Vector(cg,p5);
  v2 = Vector(cg,p6);
  area = (v1 % v2)*0.5;
  ff = (fcg+f5+f6)/3.0;
  cx += ff*area[0];
  cy += ff*area[1];
  cz += ff*area[2];
  v1 = Vector(cg,p6);
  v2 = Vector(cg,p7);
  area = (v1 % v2)*0.5;
  ff = (fcg+f6+f7)/3.0;
  cx += ff*area[0];
  cy += ff*area[1];
  cz += ff*area[2];
  v1 = Vector(cg,p7);
  v2 = Vector(cg,p4);
  area = (v1 % v2)*0.5;
  ff = (fcg+f7+f4)/3.0;
  cx += ff*area[0];
  cy += ff*area[1];
  cz += ff*area[2];

  grad = Vector(cx/vol,cy/vol,cz/vol);

  return(grad);
}

double solid_angle(Vector v1, Vector v2, Vector v3)
{
  double a, b, c, s, t, ta, tb, tc, temp, E;

  v1.normalize();
  v2.normalize();
  v3.normalize();

  a = (v1 * v2);
  b = (v1 * v3);
  c = (v2 * v3);
  a = MAX(-1.0,MIN(1.0,a));
  b = MAX(-1.0,MIN(1.0,b));
  c = MAX(-1.0,MIN(1.0,c));
  a = acos(a);
  b = acos(b);
  c = acos(c);

  s = 0.5*(a+b+c);
  t = tan(0.5*s);
  ta = tan(0.5*(s-a));
  tb = tan(0.5*(s-b));
  tc = tan(0.5*(s-c));
  temp = sqrt(MAX(0.0,t*ta*tb*tc));
  E = 4.0*atan(temp);

  return(E);
}

double solid_angle(double x0, double y0, double z0,
                   double x1, double y1, double z1,
                   double x2, double y2, double z2,
                   double x3, double y3, double z3)
{
  double a, b, c, s, t, ta, tb, tc, temp, E;
  double v1[3], v2[3], v3[3];

  v1[0] = x1-x0;
  v1[1] = y1-y0;
  v1[2] = z1-z0;
  s = MAX(1.0e-20,sqrt(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2]));
  v1[0] /= s;
  v1[1] /= s;
  v1[2] /= s;
  v2[0] = x2-x0;
  v2[1] = y2-y0;
  v2[2] = z2-z0;
  s = MAX(1.0e-20,sqrt(v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2]));
  v2[0] /= s;
  v2[1] /= s;
  v2[2] /= s;
  v3[0] = x3-x0;
  v3[1] = y3-y0;
  v3[2] = z3-z0;
  s = MAX(1.0e-20,sqrt(v3[0]*v3[0]+v3[1]*v3[1]+v3[2]*v3[2]));
  v3[0] /= s;
  v3[1] /= s;
  v3[2] /= s;

  a = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
  b = v1[0]*v3[0] + v1[1]*v3[1] + v1[2]*v3[2];
  c = v2[0]*v3[0] + v2[1]*v3[1] + v2[2]*v3[2];
  a = MAX(-1.0,MIN(1.0,a));
  b = MAX(-1.0,MIN(1.0,b));
  c = MAX(-1.0,MIN(1.0,c));
  a = acos(a);
  b = acos(b);
  c = acos(c);

  s = 0.5*(a+b+c);
  t = tan(0.5*s);
  ta = tan(0.5*(s-a));
  tb = tan(0.5*(s-b));
  tc = tan(0.5*(s-c));
  temp = sqrt(MAX(0.0,t*ta*tb*tc));
  E = 4.0*atan(temp);

  return(E);
}

void histogram(int nh, double hist[])
{
  const int nrows = 20;
  const int ncols = 80;
  char h[nrows][ncols+1];
  int bucket[ncols];
  int i, j;

  // initialize character array
  for (i=0; i < nrows; i++)
  {
    for (j=0; j < ncols; j++)
      h[i][j] = ' ';
    h[i][ncols] = '\0';
  }

  // compute ranges for buckets
  double hmin, hmax, dh;
  hmin = 1.0e20;
  hmax = -1.0e20;
  for (i=0; i < nh; i++)
  {
    hmin = MIN(hmin,hist[i]);
    hmax = MAX(hmax,hist[i]);
  }
  dh = MAX(1.0e-15,(hmax-hmin));

  // fill buckets
  for (j=0; j < ncols; j++)
    bucket[j] = 0;

  for (i=0; i < nh; i++)
  {
    j = (int)((hist[i]-hmin)/dh*ncols);
    j = MAX(0,MIN(j,ncols-1));
    bucket[j]++;
  }

  // transfer bucket count to character array
  int nmin, nmax, dn;
  nmin = bucket[0];
  nmax = bucket[0];
  for (j=1; j < ncols; j++)
  {
    nmin = MIN(nmin,bucket[j]);
    nmax = MAX(nmax,bucket[j]);
  }
  dn = MAX(1,nmax-nmin);

  for (j=0; j < ncols; j++)
  {
    i = (int)((double)(bucket[j]-nmin)/dn*nrows);
    i = MIN(nrows-1,i);
    h[i][j] = '*';
  }

  // print out histogram
  fprintf(out_f,"\nHistogram: Min function value = %g, Max function value = %g",hmin,hmax);
  fprintf(out_f,"\n%d = Max bin count",nmax);
  for (i=nrows-1; i >= 0; i--)
    fprintf(out_f,"\n|%s",h[i]);
  fprintf(out_f,"\n|");
  for (j=0; j < ncols; j++)
    fprintf(out_f,"-");
  fprintf(out_f,"\n%d = Min bin count\n",nmin);

}

int compute_Riemannian_metric(Point p1, Point p2, Point p3, Point p4, double RT[3][3])
{
  int i, j;
  double **a, *rhs, **v, *w, *x;
  double dx, dy, dz, mag;
  a = new double*[7];
  v = new double*[7];
  for (i=0; i < 7; i++)
  {
    a[i] = new double[7];
    v[i] = new double[7];
  }
  rhs = new double[7];
  w = new double[7];
  x = new double[7];

  dx = p2[0]-p1[0];
  dy = p2[1]-p1[1];
  dz = p2[2]-p1[2];
  a[1][1] = dx*dx;
  a[1][2] = 2*dx*dy;
  a[1][3] = 2*dx*dz;
  a[1][4] = dy*dy;
  a[1][5] = 2*dy*dz;
  a[1][6] = dz*dz;
  dx = p3[0]-p1[0];
  dy = p3[1]-p1[1];
  dz = p3[2]-p1[2];
  a[2][1] = dx*dx;
  a[2][2] = 2*dx*dy;
  a[2][3] = 2*dx*dz;
  a[2][4] = dy*dy;
  a[2][5] = 2*dy*dz;
  a[2][6] = dz*dz;
  dx = p4[0]-p1[0];
  dy = p4[1]-p1[1];
  dz = p4[2]-p1[2];
  a[3][1] = dx*dx;
  a[3][2] = 2*dx*dy;
  a[3][3] = 2*dx*dz;
  a[3][4] = dy*dy;
  a[3][5] = 2*dy*dz;
  a[3][6] = dz*dz;
  dx = p3[0]-p2[0];
  dy = p3[1]-p2[1];
  dz = p3[2]-p2[2];
  a[4][1] = dx*dx;
  a[4][2] = 2*dx*dy;
  a[4][3] = 2*dx*dz;
  a[4][4] = dy*dy;
  a[4][5] = 2*dy*dz;
  a[4][6] = dz*dz;
  dx = p4[0]-p2[0];
  dy = p4[1]-p2[1];
  dz = p4[2]-p2[2];
  a[5][1] = dx*dx;
  a[5][2] = 2*dx*dy;
  a[5][3] = 2*dx*dz;
  a[5][4] = dy*dy;
  a[5][5] = 2*dy*dz;
  a[5][6] = dz*dz;
  dx = p4[0]-p3[0];
  dy = p4[1]-p3[1];
  dz = p4[2]-p3[2];
  a[6][1] = dx*dx;
  a[6][2] = 2*dx*dy;
  a[6][3] = 2*dx*dz;
  a[6][4] = dy*dy;
  a[6][5] = 2*dy*dz;
  a[6][6] = dz*dz;

  rhs[1] = rhs[2] = rhs[3] = rhs[4] = rhs[5] = rhs[6] = 1.0;

  int ne = 6;

  svdcmp(a, ne, ne, w, v);

  mag = 0.0;
  for (j=1; j <= ne; j++)
    mag = MAX(mag,w[j]);
  mag *= 1.0e-06;
  for (j=1; j <= ne; j++)
  if (w[j] < mag)
    w[j] = 0.0;

  svbksb(a,w,v,ne,ne,rhs,x);

  RT[0][0] = x[1];
  RT[0][1] = RT[1][0] = x[2];
  RT[0][2] = RT[2][0] = x[3];
  RT[1][1] = x[4];
  RT[1][2] = RT[2][1] = x[5];
  RT[2][2] = x[6];
  
  for (i=0; i < 7; i++)
  {
    delete[] a[i];
    delete[] v[i];
  }
  delete[] a;
  delete[] v;
  delete[] w;
  delete[] x;
  delete[] rhs;

  return(0);
}
