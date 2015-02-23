#include <stdio.h>

#ifndef Point2D_h
#define Point2D_h
class Point2D
{
  public:
  inline Point2D(double x=0.0, double y=0.0)
  { pos[0] = x; pos[1] = y; }
  ~Point2D() { }
  void print(FILE *outf)
  {
       fprintf(outf,"\nPoint (x,y)= (%.12g,%.12g)",pos[0],pos[1]);
  }
  inline Point2D operator + (const Point2D &) const;
  inline Point2D operator - (const Point2D &) const;
  inline Point2D &operator += (const Point2D &);
  inline Point2D &operator -= (const Point2D &);
  inline Point2D &operator += (double);
  inline Point2D &operator -= (double);
  inline Point2D &operator *= (double);
  inline Point2D &operator /= (double);
  inline Point2D operator * (double) const;
  inline Point2D operator / (double) const;
  inline double & operator () (int);
  inline double operator () (int) const;
  inline double & operator [] (int);

  double pos[2];
};

inline double & Point2D::operator () (int i) 
{ return pos[i]; }

inline double Point2D::operator () (int i) const
{ return pos[i]; }

inline double & Point2D::operator [] (int i) 
{ return pos[i]; }

inline Point2D Point2D::operator * (double other) const 
{ return Point2D(pos[0]*other, pos[1]*other); }

inline Point2D Point2D::operator / (double other) const 
{
#ifdef _DEBUG_
    if (fabs(other) <=1e-20) throw new mException (__LINE__,__FILE__);
#endif
    return Point2D(pos[0]/other, pos[1]/other);
}

inline Point2D Point2D::operator + (const Point2D &other) const 
{ return Point2D(pos[0]+other.pos[0], pos[1]+other.pos[1]); }

inline Point2D Point2D::operator - (const Point2D &other) const 
{ return Point2D(pos[0]-other.pos[0], pos[1]-other.pos[1]); }

inline Point2D& Point2D::operator += (const Point2D &other) 
{
	for (int i = 0; i < 2; i++)
          {
            pos[i]+=other.pos[i];
          }
	return *this;
}

inline Point2D& Point2D::operator -= (const Point2D &other) 
{
	for (int i = 0; i < 2; i++)
          {
            pos[i]-=other.pos[i];
          }
	return *this;
}

inline Point2D& Point2D::operator += (double other) 
{
	for (int i = 0; i < 2; i++)
          {
            pos[i]+=other;
          }
	return *this;
}

inline Point2D& Point2D::operator -= (double other) 
{
	for (int i = 0; i < 2; i++)
          {
            pos[i]-=other;
          }
	return *this;
}

inline Point2D& Point2D::operator *= (double other) 
{
	for (int i = 0; i < 2; i++)
          {
            pos[i]*=other;
          }
	return *this;
}

inline Point2D& Point2D::operator /= (double other) 
{
#ifdef _DEBUG_
	if (fabs(other) <=1e-20) throw new mException (__LINE__,__FILE__);
#endif
	for (int i = 0; i < 2; i++)
          {
            pos[i]/=other;
          }
	return *this;
}

#endif
