#include <stdio.h>
#include <math.h>
#include "space.h"
#include "Point.h"

#ifndef Vector_h
#define Vector_h
class Vector
{
 double vec[SPACE];
 public:
  #if SPACE == 2
  inline Vector(double dx=0.0, double dy=0.0)
  { vec[0] = dx; vec[1] = dy; }
  #else
  inline Vector(double dx=0.0, double dy=0.0, double dz=0.0)
  { vec[0] = dx; vec[1] = dy; vec[2] = dz; }
  #endif
  inline Vector(const Point &from, const Point &to);
  #if SPACE == 2
  void print(FILE *fp)
  { fprintf(fp,"\nVector (x,y)= (%.12g,%.12g)",vec[0],vec[1]); }
  #else
  void print(FILE *fp)
  { fprintf(fp,"\nVector (x,y,z)= (%.12g,%.12g,%.12g)",vec[0],vec[1],vec[2]); }
  #endif
  inline Vector operator + ( const Vector &v);
  inline Vector operator - ( const Vector &v);
  inline Vector &operator += ( const Vector &v);
  inline Vector &operator -= ( const Vector &v);
  inline Vector &operator += ( const double &s);
  inline Vector &operator -= ( const double &s);
  inline Vector &operator *= ( const double &s);
  inline Vector &operator /= ( const double &s);
  #if SPACE == 2
  inline double operator * ( const Vector &v)
  { return double(vec[0]*v.vec[0]+vec[1]*v.vec[1]); }
  #else
  inline double operator * ( const Vector &v)
  { return double(vec[0]*v.vec[0]+vec[1]*v.vec[1]+vec[2]*v.vec[2]); }
  #endif
  #if SPACE == 2
  inline double operator % ( const Vector &v);
  #else
  inline Vector operator % ( const Vector &v);
  #endif
  inline Vector operator * ( const double &s);
  inline Vector operator / ( const double &s);
  inline double operator () (int i) const;
  inline double &operator () (int i);
  inline double &operator [] (int i);
  inline double magnitude();
  inline void normalize();
  ~Vector(){};

};

inline Vector::Vector(const Point &from, const Point &to)
{
  for (int i = 0; i < SPACE; i++)
          {
            vec[i] = to(i) - from(i);
          }
}
inline double &Vector::operator () (int i)
{ return vec[i]; }
inline double Vector::operator () (int i) const
{ return vec[i]; }
inline double &Vector::operator [] (int i)
{ return vec[i]; }
#if SPACE == 2
inline double Vector::operator % ( const Vector &v)
{
  return (vec[0]*v.vec[1]-vec[1]*v.vec[0]);
}
#else
inline Vector Vector::operator % ( const Vector &v)
{
  return Vector(vec[1]*v.vec[2]-vec[2]*v.vec[1],
                vec[2]*v.vec[0]-vec[0]*v.vec[2],
                vec[0]*v.vec[1]-vec[1]*v.vec[0]);
}
#endif
#if SPACE == 2
inline Vector Vector::operator +( const Vector &v)
{ return Vector(vec[0]+v.vec[0], vec[1]+v.vec[1]); }
inline Vector Vector::operator -( const Vector &v)
{ return Vector(vec[0]-v.vec[0], vec[1]-v.vec[1]); }
inline Vector Vector::operator *( const double &s)
{ return Vector(s*vec[0], s*vec[1]); }
inline Vector Vector::operator /( const double &s)
{ return Vector(vec[0]/s, vec[1]/s); }
inline double Vector::magnitude()
{ return double(sqrt(vec[0]*vec[0] + vec[1]*vec[1])); }
#else
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
#endif
inline Vector &Vector::operator += (const Vector &v)
{
  for (int i = 0; i < SPACE; i++)
          {
            vec[i]+=v.vec[i];
          }
  return *this;
}
inline Vector &Vector::operator -= (const Vector &v)
{
  for (int i = 0; i < SPACE; i++)
          {
            vec[i]-=v.vec[i];
          }
  return *this;
}
inline Vector &Vector::operator += (const double &s)
{
  for (int i = 0; i < SPACE; i++)
          {
            vec[i]+=s;
          }
  return *this;
}
inline Vector &Vector::operator -= (const double &s)
{
  for (int i = 0; i < SPACE; i++)
          {
            vec[i]-=s;
          }
  return *this;
}
inline Vector &Vector::operator *= (const double &s)
{
  for (int i = 0; i < SPACE; i++)
          {
            vec[i]*=s;
          }
  return *this;
}
inline Vector &Vector::operator /= (const double &s)
{
  for (int i = 0; i < SPACE; i++)
          {
            vec[i]/=s;
          }
  return *this;
}

// we do not need to worry if mag less than machine eps
// if so, vector is 0 anyways, and if we divide nearly 0 by zero,
// we get a huge number and not a normalized vector like we want 
#if SPACE == 2
inline void Vector::normalize()
{
  double mag = sqrt(vec[0]*vec[0]+vec[1]*vec[1]);
  if (mag > 1.0e-20) { vec[0] = vec[0] / mag; vec[1] = vec[1] / mag; }
}
#else
inline void Vector::normalize()
{
  double mag = sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
  if (mag > 1.0e-20) { vec[0] = vec[0] / mag; vec[1] = vec[1] / mag; vec[2] = vec[2] / mag; }
}
#endif

#if SPACE == 3
inline double scalar_triple_product(Vector a, Vector b, Vector c)
{
  double d;
  d = a[0]*b[1]*c[2]+a[1]*b[2]*c[0]+a[2]*b[0]*c[1]-
      a[2]*b[1]*c[0]-a[1]*b[0]*c[2]-a[0]*b[2]*c[1];
  return(d);
}
#endif

#endif
