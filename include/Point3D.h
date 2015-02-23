#include <stdio.h>

#ifndef Point_h
#define Point_h
class Point
{
  public:
  inline Point(double x=0.0, double y=0.0, double z=0.0)
  { pos[0] = x; pos[1] = y; pos[2] = z; }
  virtual ~Point() { }
  virtual void print(FILE *fp)
  { fprintf(fp,"\nPoint (x,y,z)= (%g,%g,%g)",pos[0],pos[1],pos[2]); }
  inline Point operator + (const Point &) const;
  inline Point operator - (const Point &) const;
  inline Point &operator += (const Point &);
  inline Point &operator -= (const Point &);
  inline Point &operator += (double);
  inline Point &operator -= (double);
  inline Point operator * (double) const;
  inline Point operator / (double) const;
  inline Point& operator *= (double);
  inline Point& operator /= (double);
  inline double & operator () (int);
  inline double operator () (int) const;
  inline double & operator [] (int);

  double pos[3];
};

inline double & Point::operator () (int i) 
{ return pos[i]; }

inline double Point::operator () (int i) const
{ return pos[i]; }

inline double & Point::operator [] (int i) 
{ return pos[i]; }

inline Point Point::operator * (double other) const 
{ return Point(pos[0]*other, pos[1]*other, pos[2]*other); }

inline Point Point::operator / (double other) const 
{
#ifdef _DEBUG_
	if (fabs(other) <=1e-20) throw new mException (__LINE__,__FILE__);
#endif
  return Point(pos[0]/other, pos[1]/other,pos[2]/other);
}

inline Point Point::operator + (const Point &other) const 
{ return Point(pos[0]+other.pos[0], pos[1]+other.pos[1], pos[2]+other.pos[2]); }

inline Point Point::operator - (const Point &other) const 
{ return Point(pos[0]-other.pos[0], pos[1]-other.pos[1], pos[2]-other.pos[2]); }

inline Point& Point::operator += (const Point &other) 
{
	pos[0]+=other.pos[0]; 
	pos[1]+=other.pos[1];
	pos[2]+=other.pos[2];
	return *this;
}

inline Point& Point::operator -= (const Point &other) 
{
	pos[0]-=other.pos[0]; 
	pos[1]-=other.pos[1];
	pos[2]-=other.pos[2];
	return *this;
}

inline Point& Point::operator += (double other) 
{
	pos[0]+=other; 
	pos[1]+=other;
	pos[2]+=other;
	return *this;
}

inline Point& Point::operator -= (double other) 
{
	pos[0]-=other; 
	pos[1]-=other;
	pos[2]-=other;
	return *this;
}

inline Point& Point::operator *= (double other) 
{
	pos[0]*=other; 
	pos[1]*=other;
	pos[2]*=other;
	return *this;
}

inline Point& Point::operator /= (double other) 
{
#ifdef _DEBUG_
	if (fabs(other) <=1e-20) throw new mException (__LINE__,__FILE__);
#endif
	pos[0]/=other; 
	pos[1]/=other;
	pos[2]/=other;
	return *this;
}

#endif
