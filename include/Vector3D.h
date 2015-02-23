#include <stdio.h>
#include <math.h>
#include "Point.h"

#ifndef Vector_h
#define Vector_h
class Vector
{
 double vec[3];
 public:
  inline Vector(double dx=0.0, double dy=0.0, double dz=0.0)
  { vec[0] = dx; vec[1] = dy; vec[2] = dz; }
  inline Vector(const Point &from, const Point &to);
  void print(FILE *fp)
  { fprintf(fp,"\nVector (x,y,z)= (%g,%g,%g)",vec[0],vec[1],vec[2]); }
  inline Vector operator + ( const Vector &v);
  inline Vector operator - ( const Vector &v);
  inline Vector &operator += ( const Vector &v);
  inline Vector &operator -= ( const Vector &v);
  inline Vector &operator += ( const double &s);
  inline Vector &operator -= ( const double &s);
  inline Vector &operator *= ( const double &s);
  inline Vector &operator /= ( const double &s);
  inline double operator * ( const Vector &v)
  { return double(vec[0]*v.vec[0]+vec[1]*v.vec[1]+vec[2]*v.vec[2]); }
  inline Vector operator % ( const Vector &v);
  inline Vector operator * ( const double &s);
  inline Vector operator / ( const double &s);
  inline double operator () (int i) const;
  inline double &operator () (int i);
  inline double &operator [] (int i);
  inline double magnitude();
  inline void normalize();
  inline ~Vector(){};

};

inline Vector::Vector(const Point &from, const Point &to)
{
  vec[0] = to(0) - from(0);
  vec[1] = to(1) - from(1);
  vec[2] = to(2) - from(2);
}
inline double &Vector::operator () (int i)
{ return vec[i]; }
inline double Vector::operator () (int i) const
{ return vec[i]; }
inline double &Vector::operator [] (int i)
{ return vec[i]; }
inline Vector Vector::operator % ( const Vector &v)
{
  return Vector(vec[1]*v.vec[2]-vec[2]*v.vec[1],
                vec[2]*v.vec[0]-vec[0]*v.vec[2],
                vec[0]*v.vec[1]-vec[1]*v.vec[0]);
}
inline Vector Vector::operator +( const Vector &v)
{ return Vector(vec[0]+v.vec[0], vec[1]+v.vec[1], vec[2]+v.vec[2]); }
inline Vector Vector::operator -( const Vector &v)
{ return Vector(vec[0]-v.vec[0], vec[1]-v.vec[1], vec[2]-v.vec[2]); }
inline Vector Vector::operator *( const double &s)
{ return Vector(s*vec[0], s*vec[1], s*vec[2]); }
inline Vector Vector::operator /( const double &s)
{ return Vector(vec[0]/s, vec[1]/s, vec[2]/s); }
inline double Vector::magnitude()
{ return double(sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2])); }
inline Vector &Vector::operator += (const Vector &v)
{
  vec[0]+=v.vec[0];
  vec[1]+=v.vec[1];
  vec[2]+=v.vec[2];
  return *this;
}
inline Vector &Vector::operator -= (const Vector &v)
{
  vec[0]-=v.vec[0];
  vec[1]-=v.vec[1];
  vec[2]-=v.vec[2];
  return *this;
}
inline Vector &Vector::operator += (const double &s)
{
  vec[0]+=s;
  vec[1]+=s;
  vec[2]+=s;
  return *this;
}
inline Vector &Vector::operator -= (const double &s)
{
  vec[0]-=s;
  vec[1]-=s;
  vec[2]-=s;
  return *this;
}
inline Vector &Vector::operator *= (const double &s)
{
  vec[0]*=s;
  vec[1]*=s;
  vec[2]*=s;
  return *this;
}
inline Vector &Vector::operator /= (const double &s)
{
  vec[0]/=s;
  vec[1]/=s;
  vec[2]/=s;
  return *this;
}
inline void Vector::normalize()
{
  double mag = sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
  if (mag > 1.0e-20) { vec[0] = vec[0] / mag; vec[1] = vec[1] / mag; vec[2] = vec[2] / mag; }
}

inline double scalar_triple_product(Vector a, Vector b, Vector c)
{
  double d;
  d = a[0]*b[1]*c[2]+a[1]*b[2]*c[0]+a[2]*b[0]*c[1]-
      a[2]*b[1]*c[0]-a[1]*b[0]*c[2]-a[0]*b[2]*c[1];
  return(d);
}

#endif
