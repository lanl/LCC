#include<fstream>
#include<string>
#include<iostream>
#include<iomanip>
#include<cmath>
#include<vector>
#include<random>

using namespace std;


void grabber(string input, string &first, string &last, char mid){
  string hold;
  hold = input;
  int index = 0;
  bool midthere = false;
  for(int i = 0 ; i < input.length(); i++){
    if(hold[i] == mid){
      index = i;
      midthere = true;
    }
  }
  if(midthere == false)
    first = input;
  else{
    last = input.substr(index+1);
    first = input.substr(0, index);
  }
}

void setter(string first, string last, double& nx, double& ny, double& nz, bool& random, string& lattice, int& t, double& tmax, double& tstart, double& tol, double& tstep) {
  if(first == "random"){
    if(last == "n" || last == "false" || last == "no" || last == "off")
      random = false;
    else
      random = true;
  }
  if(first == "nx")
    nx = stod(last);
  if(first == "ny")
    ny = stod(last);
  if(first == "nz")
    nz = stod(last);
  if(first == "lattice")
    lattice = last;
  if(first == "t")
    t = stoi(last);
  if(first == "tmax")
    tmax = stod(last);
  if(first == "tstart"){
    if(stod(last) == 0)
      tstart = 1;
    else
      tstart = stod(last);
  }
    if(first == "tol")
    tol = stod(last);
    if(first == "tstep")
    tstep = stod(last);
}

class Vector
{
    public:
  double x;
  double y;
  double z;

        // Default constructor: create a vector whose
        // x, y, z components are all zero.
 Vector()
    : x(0.0)
    , y(0.0)
    , z(0.0)
  {
  }

  // This constructor initializes a vector to any desired component values. 
  Vector(double _x, double _y, double _z)
    : x(_x)
    , y(_y)
    , z(_z)
  {
  }
  // Returns the square of the magnitude of this vector.
  // This is more efficient than computing the magnitude itself,
  // and is just as good for comparing two vectors to see which
  // is longer or shorter.                                                                                                                                                                                  
  const double MagnitudeSquared() const
  {
    return (x*x) + (y*y) + (z*z);
  }

  const double Magnitude() const
  {
    return sqrt(MagnitudeSquared());
  }

  const Vector UnitVector() const
  {
    const double mag = Magnitude();
    return Vector(x/mag, y/mag, z/mag);
}

  Vector& operator *= (const double factor)
  {
    x *= factor;
    y *= factor;
    z *= factor;
    return *this;
  }

  Vector& operator += (const Vector& other)
  {
    x += other.x;
    y += other.y;
    z += other.z;
    return *this;
  }
};
//-------------------------------------------------------------
inline Vector operator + (const Vector &a, const Vector &b){
  return Vector(a.x + b.x, a.y + b.y, a.z + b.z);
}

inline Vector operator - (const Vector &a, const Vector &b){
  return Vector(a.x - b.x, a.y - b.y, a.z - b.z);
}

inline Vector operator - (const Vector& a){
  return Vector(-a.x, -a.y, -a.z);
}

inline double DotProduct (const Vector& a, const Vector& b){
  return (a.x*b.x) + (a.y*b.y) + (a.z*b.z);
}

inline Vector CrossProduct (const Vector& a, const Vector& b){
  return Vector(
                (a.y * b.z) - (a.z * b.y),
                (a.z * b.x) - (a.x * b.z),
                (a.x * b.y) - (a.y * b.x));
}

inline Vector operator * (double s, const Vector& v){
  return Vector(s*v.x, s*v.y, s*v.z);
}

inline Vector operator / (const Vector& v, double s){
  return Vector(v.x/s, v.y/s, v.z/s);
}
//-------------------------------------------------------------



