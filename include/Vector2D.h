#include <stdio.h>
#include <math.h>
#include "Point2D.h"

#ifndef Vector2D_h
#define Vector2D_h
class Vector2D
{
 double vec[2];
 public:
  inline Vector2D(double dx=0.0, double dy=0.0)
  { vec[0] = dx; vec[1] = dy; }
  inline Vector2D(const Point2D &from, const Point2D &to);
  void print(FILE *outf)
  { fprintf(outf,"\nVector (x,y)= (%.12g,%.12g)",vec[0],vec[1]); }
  inline Vector2D operator + ( const Vector2D &v);
  inline Vector2D operator - ( const Vector2D &v);
  inline Vector2D &operator += ( const Vector2D &v);
  inline Vector2D &operator -= ( const Vector2D &v);
  inline Vector2D &operator += ( const double &s);
  inline Vector2D &operator -= ( const double &s);
  inline Vector2D &operator *= ( const double &s);
  inline Vector2D &operator /= ( const double &s);
  inline double operator * ( const Vector2D &v)
  { return double(vec[0]*v.vec[0]+vec[1]*v.vec[1]); }
  inline double operator % ( const Vector2D &v);
  inline Vector2D operator * ( const double &s);
  inline Vector2D operator / ( const double &s);
  inline double operator () (int i) const;
  inline double &operator () (int i);
  inline double &operator [] (int i);
  inline double magnitude();
  inline void normalize();
  ~Vector2D(){};

};

inline Vector2D::Vector2D(const Point2D &from, const Point2D &to)
{
  for (int i = 0; i < 2; i++)
          {
            vec[i] = to(i) - from(i);
          }
}
inline double &Vector2D::operator () (int i)
{ return vec[i]; }
inline double Vector2D::operator () (int i) const
{ return vec[i]; }
inline double &Vector2D::operator [] (int i)
{ return vec[i]; }
inline double Vector2D::operator % ( const Vector2D &v)
{
  return (vec[0]*v.vec[1]-vec[1]*v.vec[0]);
}
inline Vector2D Vector2D::operator +( const Vector2D &v)
{ return Vector2D(vec[0]+v.vec[0], vec[1]+v.vec[1]); }
inline Vector2D Vector2D::operator -( const Vector2D &v)
{ return Vector2D(vec[0]-v.vec[0], vec[1]-v.vec[1]); }
inline Vector2D Vector2D::operator *( const double &s)
{ return Vector2D(s*vec[0], s*vec[1]); }
inline Vector2D Vector2D::operator /( const double &s)
{ return Vector2D(vec[0]/s, vec[1]/s); }
inline double Vector2D::magnitude()
{ return double(sqrt(vec[0]*vec[0] + vec[1]*vec[1])); }
inline Vector2D &Vector2D::operator += (const Vector2D &v)
{
  for (int i = 0; i < 2; i++)
          {
            vec[i]+=v.vec[i];
          }
  return *this;
}
inline Vector2D &Vector2D::operator -= (const Vector2D &v)
{
  for (int i = 0; i < 2; i++)
          {
            vec[i]-=v.vec[i];
          }
  return *this;
}
inline Vector2D &Vector2D::operator += (const double &s)
{
  for (int i = 0; i < 2; i++)
          {
            vec[i]+=s;
          }
  return *this;
}
inline Vector2D &Vector2D::operator -= (const double &s)
{
  for (int i = 0; i < 2; i++)
          {
            vec[i]-=s;
          }
  return *this;
}
inline Vector2D &Vector2D::operator *= (const double &s)
{
  for (int i = 0; i < 2; i++)
          {
            vec[i]*=s;
          }
  return *this;
}
inline Vector2D &Vector2D::operator /= (const double &s)
{
  for (int i = 0; i < 2; i++)
          {
            vec[i]/=s;
          }
  return *this;
}

// we do not need to worry if mag less than machine eps
// if so, vector is 0 anyways, and if we divide nearly 0 by zero,
// we get a huge number and not a normalized vector like we want 
inline void Vector2D::normalize()
{
  double mag = sqrt(vec[0]*vec[0]+vec[1]*vec[1]);
  if (mag > 1.0e-20) { vec[0] = vec[0] / mag; vec[1] = vec[1] / mag; }
}

#endif
