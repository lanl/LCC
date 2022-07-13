//Andrew Alvarado script for generating encapsulating orthogonal vectors
//uses functions in func_v.h
#include"func_v.h"
#include<ctime>

using namespace std;
int main(int argc, char** argv){
  ifstream infile;
  ofstream outfile, outfile2;
  string input, first, last, temp, temp1;
  string lattice = "bcc";
  int count = 0;
  string filename = "";
  bool random = false;
  double nx=1, ny=1, nz=1;
  int t = 1;
  double tmax = 2;
  double tstart = 1;
  double tstep = 0.01;
  double tol = 0.0001;
  if(argc <= 1){
    cout << "start vector n={nx,ny,xz} needed, exp: nx= , ny= , nz= | or random=y, optional t={0, 1, 2,...} , lattice={fcc,bcc} (no brackets), tmax={2.0}, tstart={0.1}, tstep={0.01}, tol={0.0001}" << endl;
      return 0;
  }

  default_random_engine generator(random_device{}());
  mt19937 gen;
  uniform_real_distribution<double> unifrnd(0.0, 1.0);
  
  if(argc > 1){
    for(int i = 1; i < argc; i++){
      input = argv[i];
      grabber(input, first, last, '=');
      setter(first, last, nx, ny, nz, random, lattice, t, tmax, tstart, tol, tstep);
    }
  }
  //----------------initialization--------------
  Vector n1;
  if(random == true)
    n1 = Vector(unifrnd(generator), unifrnd(generator), unifrnd(generator));
  else
    n1 = Vector(nx,ny,nz);
  Vector n2, n3;
  //------------find orthogonal vector----------
  if(n1.x == 0){ //n1.x is zero
    if(n1.y == 0){
      if(n1.z == 0){ //0 vector
	cout << "start vector is {0,0,0}" << endl;
	return 0;
      } //end of 0 vector
      else{ //n1.x=0 & n1.y=0
	n2.x = 0.0;
	n2.y = n1.z;
	n2.z = 0.0;
      }
    } //end of n1.x=0 & n1.y=0
    else{ //n1.y non-zero case
      if(n1.z == 0){ //x & z are zero, y is non-zero;
	n2.x = n1.y;
	n2.y = 0.0;
	n2.z = 0.0;
      } //end x & z are zero, y is non-zero;
      else{// x is zero, y is non-zero, z is non-zero
	n2.x = 0.0;
	n2.y = -n1.z/n1.y;
	n2.z = 1.0;
      } // end x is zero, y is non-zero, z is non zero
    } //end n1.y is non-zero
  }
  else{ //n1.x is non-zero
    if(n1.y == 0){
      if(n1.z == 0){ // n1.y=0 & n1.z=0
	n2.x = 0.0;
	n2.y = 1.0;
	n2.z = 0.0;
      } // end n1.y=0 & n1.z=0
      else{ // n1.y is zero
	n2.x = 1.0;
	n2.y = 0.0; 
	n2.z = -n1.x/n1.z;
      } // end n1.y is zero
    }
    else{ // n1.y is non-zero
      if(n1.z == 0){ // n1.z is zero
	n2.x = -n1.y/n1.x;
	n2.y = 1.0;
	n2.z = 0.0;
      } // end n1.z is zero
      else{ // all non-zero
	n2.x = -n1.y/n1.x;
	n2.y = 1;
	n2.z = 0;
      }
    }
  }  
  //---------find orthogonal vector end---------
  //n3 = CrossProduct(n1, n2);  //or
  if(n1.x == 0 && n1.y ==0){
    n3.x = 1.0;
    n3.y = 0;
    n3.z = 0;
  }
  else if(n1.x == 0 && n1.z == 0){
    n3.x = 0;
    n3.y = 0;
    n3.z = 1.0;
  }
 else if(n1.y == 0 && n1.z == 0){
    n3.x= 0;
    n3.y = 0;
    n3.z = 1.0;
  }
  else if(n1.x == 0 && n1.y != 0 && n1.z !=0){
    n3.x = 1.0;
    n3.y = ((n1.z*n1.z*n1.x)/(n1.y*n1.y*n1.y))/( ((n1.z*n1.z)/(n1.y*n1.y)) +1.0 ) - n1.x/n1.y;
    n3.z = ((-n1.z*n1.x)/(n1.y*n1.y))/(n1.z*n1.z/(n1.y*n1.y)+1);  
  }
  else if(n1.x != 0 && n1.y == 0 && n1.z !=0 ){
    n3.x = ((-n1.x*n1.y)/(n1.z*n1.z))/(n1.x*n1.x/(n1.z*n1.z)+1);
    n3.y = 1.0; 
    n3.z = ((n1.x*n1.x*n1.y)/(n1.z*n1.z*n1.z))/( ((n1.x*n1.x)/(n1.z*n1.z)) +1.0 ) - n1.y/n1.z;

  }
  else{
  n3.x = ((n1.y*n1.y*n1.z)/(n1.x*n1.x*n1.x))/( ((n1.y*n1.y)/(n1.x*n1.x)) +1.0 ) - n1.z/n1.x;
  n3.y = ((-n1.y*n1.z)/(n1.x*n1.x))/(n1.y*n1.y/(n1.x*n1.x)+1);  
  n3.z = 1.0; //fails for n1=(1, 0, 1)
  }
  //------------end initialization---------------
  cout << "n0: " << setw(3) << n1.x << setw(12) << n1.y << setw(12) << n1.z << endl;
  cout << "n1: " << setw(3) << n2.x << setw(12) << n2.y << setw(12) << n2.z << endl;
  cout << "n2: " << setw(3) << n3.x << setw(12) << n3.y << setw(12) << n3.z << endl; 
  cout << "n1•n2 = " << DotProduct(n1,n2) << endl;
  cout << "n1•n3 = " << DotProduct(n1,n3) << endl;
  cout << "n2•n3 = " << DotProduct(n2,n3) << endl;
  
  cout << fixed << right << setprecision(5);
  outfile.open("output.txt");
  outfile << "#t1 \t eps \t  a \t  b \t c" << endl;
  outfile2.open("output2.txt");
  outfile2 << "#angle eps" << endl;
  //find t such that error is lowest for given lattice vectors
  Vector a1, a2, a3, m1, m2, m3;//a's are lattice vectors m's are adjusted for matrix
  double a, b, c;
  double detM, eps;
  if(lattice=="fcc"){
    a1 = t*Vector(0.0, 0.5, 0.5);
    a2 = t*Vector(0.5, 0.0, 0.5);
    a3 = t*Vector(0.5, 0.5, 0.0);
  }
  else if (lattice == "bcc"){
    a1 = t*Vector(1.0, 0.0, 0.0);
    a2 = t*Vector(0.0, 1.0, 0.0);
    a3 = t*Vector(0.5, 0.5, 0.5);
  }
  else if (lattice == "cubic"){
    a1 = t*Vector(1.0, 0.0, 0.0);
    a2 = t*Vector(0.0, 1.0, 0.0);
    a3 = t*Vector(0.0, 0.0, 1.0);
  }    
  else if (lattice == "read"){
    infile.open("latticefile.dat");
    if(infile.good() != true){
      cout << "latticefile.dat" << " not found" << endl;
      return 0;
    }
    infile >> a1.x >> a1.y >> a1.z;
    infile >> a2.x >> a2.y >> a2.z;
    infile >> a3.x >> a3.y >> a3.z;
    infile.close();
  }
  m1 = Vector(a1.x, a2.x, a3.x);
  m2 = Vector(a1.y, a2.y, a3.y);
  m3 = Vector(a1.z, a2.z, a3.z);
  
  cout << endl << "Lattice in use: " << endl;
  cout << setw(12) << m1.x << setw(12) << m1.y << setw(12) << m1.z << endl;
  cout << setw(12) << m2.x << setw(12) << m2.y << setw(12) << m2.z << endl;
  cout << setw(12) << m3.x << setw(12) << m3.y << setw(12) << m3.z << endl << endl;
  
  double p, q, r, loweps = 1000, lowalpha, low[3], lowt[3]; //storing lowest eps values;
  double t1=1, t2=1, t3=1;
  Vector RotMat[3];
  Vector tempvec1, tempvec2, tempvec3;
  Vector lowest1, lowest2, lowest3;
  
  p = n1.UnitVector().x; //set rotational axis equivalent to n1.
  q = n1.UnitVector().y; 
  r = n1.UnitVector().z;
  
  //inverse matrix M-1 = 1/det(M) adj(intermediate_M)
  detM = (m1.x*(m2.y*m3.z-m3.y*m2.z) - m2.x*(m1.y*m3.z-m3.y*m1.z) + m3.x*(m1.y*m2.z-m2.y*m1.z)); 
  if(detM == 0){
    cout << "not invertible" << endl;
    return 0;
  }
  // |   a2.y*a3.z-a3.y*a2.z  -(a1.y*a3.z-a3.y*a1.z)   a1.y*a2.z-a2.y*a1.z  |
  // | -(a2.x*a3.z-a3.x*a2.z)   a1.x*a3.z-a3.x*a1.z  -(a1.x*a2.z-a2.x*a1.z) |
  // |   a2.x*a3.y-a3.x*a2.y  -(a1.x*a3.y-a3.x*a1.y)   a1.x*a2.y-a2.x*a1.y  | 
  /*  Vector Minv[3] = {
    Vector(a2.y*a3.z-a3.y*a2.z, -(a1.y*a3.z-a3.y*a1.z), a1.y*a2.z-a2.y*a1.z),
    Vector(-(a2.x*a3.z-a3.x*a2.z), a1.x*a3.z-a3.x*a1.z, -(a1.x*a2.z-a2.x*a1.z)),
    Vector( a2.x*a3.y-a3.x*a2.y,-(a1.x*a3.y-a3.x*a1.y) , a1.x*a2.y-a2.x*a1.y)
    };*/

  Vector Minv[3] = {
    Vector(m2.y*m3.z-m3.y*m2.z, -(m2.x*m3.z-m3.x*m2.z), m2.x*m3.y-m3.x*m2.y),
    Vector(-(m1.y*m3.z-m3.y*m1.z), m1.x*m3.z-m3.x*m1.z, -(m1.x*m3.y-m3.x*m1.y)),
    Vector( m1.y*m2.z-m2.y*m1.z,-(m1.x*m2.z-m2.x*m1.z), m1.x*m2.y-m2.x*m1.y)
  };
  
  Minv[0] = 1/detM * Minv[0];
  Minv[1] = 1/detM * Minv[1];
  Minv[2] = 1/detM * Minv[2];
  cout << "Inverse Matrix: " << endl;
  cout << setw(12) << Minv[0].x << setw(12) << Minv[0].y << setw(12) << Minv[0].z << endl;
  cout << setw(12) << Minv[1].x << setw(12) << Minv[1].y << setw(12) << Minv[1].z << endl;
  cout << setw(12) << Minv[2].x << setw(12) << Minv[2].y << setw(12) << Minv[2].z << endl << endl;
  cout << "M•M-1 : " << endl;
  cout << setw(12) << DotProduct(m1, Minv[0]) << setw(12) << DotProduct(m1, Minv[1]) << setw(12) << DotProduct(m1, Minv[2]) << endl;
  cout << setw(12) << DotProduct(m2, Minv[0]) << setw(12) << DotProduct(m2, Minv[1]) << setw(12) << DotProduct(m2, Minv[2]) << endl;
  cout << setw(12) << DotProduct(m3, Minv[0]) << setw(12) << DotProduct(m3, Minv[1]) << setw(12) << DotProduct(m3, Minv[2]) << endl << endl;
  //find length for plane normal
  //set up inverse matrix find a b and c
  //error function eps = |whole - part|           
  for(double l=tstart; l<tmax; l=l+tstep){
    tempvec1 = l*n1;
    tempvec2 = n2;
    tempvec3 = n3;
    a = DotProduct(tempvec1,Minv[0]);
    b = DotProduct(tempvec1,Minv[1]);
    c = DotProduct(tempvec1,Minv[2]);
    eps = abs(a-round(a))+abs(b-round(b))+abs(c-round(c));
    if(eps < loweps-tol){
      loweps = eps;
      lowt[0] = l;
     }
    outfile << l << " " << eps << " " <<  a << " " << b << " " << c << endl; 
  } // end t1
  outfile << endl;
  t1 = lowt[0];
  low[0] = loweps;
  loweps = 1000;
  for(double l=tstart; l<tmax; l=l+tstep){
    tempvec1 = n1;
    tempvec2 = l*n2;
    tempvec3 = n3;
    a = DotProduct(tempvec2,Minv[0]);
    b = DotProduct(tempvec2,Minv[1]);
    c = DotProduct(tempvec2,Minv[2]);
    eps = abs(a-round(a))+abs(b-round(b))+abs(c-round(c));
    if(eps < loweps-tol){
      loweps = eps;
      lowt[1] = l;
    }
    outfile << l << " " << eps << " " <<  a << " " << b << " " << c << endl;
  } //end t2
  outfile << endl;
  t2 = lowt[1];
  low[1] = loweps;
  loweps = 1000;
  for(double l=tstart; l<tmax; l=l+tstep){
    tempvec1 = n1;
    tempvec2 = n2;
    tempvec3 = l*n3;
    a = DotProduct(tempvec3,Minv[0]);
    b = DotProduct(tempvec3,Minv[1]);
    c = DotProduct(tempvec3,Minv[2]);
    eps = abs(a-round(a))+abs(b-round(b))+abs(c-round(c));
    if(eps < loweps-tol){
      loweps = eps;
      lowt[2] = l;
    }
    outfile << l << " " << eps << " " <<  a << " " << b << " " << c << endl;
  } // end t3
  outfile << endl;
  t3 = lowt[2];
  low[2] = loweps;
  loweps=1000;
  //rotation
  for(double alpha = 0; alpha <= 180*M_PI/180; (alpha=alpha+0.50*M_PI/180)){
    //rotate matrix by angle alpha around axis n1 (unitvector).
    RotMat[0] = Vector((cos(alpha)+p*p*(1-cos(alpha))),
		       (p*q*(1-cos(alpha))-r*sin(alpha)),
		       (p*r*(1-cos(alpha))+q*sin(alpha)));
    
    RotMat[1] = Vector((q*p*(1-cos(alpha))+r*sin(alpha)),
		       (cos(alpha)+q*q*(1-cos(alpha))),
		       (q*r*(1-cos(alpha))-p*sin(alpha)));
    
    RotMat[2] = Vector((r*p*(1-cos(alpha))-q*sin(alpha)),
		       (r*q*(1-cos(alpha))+p*sin(alpha)),
		       (cos(alpha)+r*r*(1-cos(alpha))));
    //rotated n2 and n3 by alpha around n1
    tempvec1 = t1*n1;
    tempvec2 = t2*Vector(DotProduct(RotMat[0], n2),
			 DotProduct(RotMat[1], n2),
			 DotProduct(RotMat[2], n2));
    tempvec3 = t3*Vector(DotProduct(RotMat[0], n3),
			 DotProduct(RotMat[1], n3),
			 DotProduct(RotMat[2], n3));
    
    a = (1/detM)*DotProduct(tempvec2,Minv[0]);
    b = (1/detM)*DotProduct(tempvec2,Minv[1]);
    c = (1/detM)*DotProduct(tempvec2,Minv[2]);
    
    //error function |whole - part|
    eps = abs(a-round(a))+abs(b-round(b))+abs(b-round(b));
    outfile2 << alpha*180/M_PI << setw(12) << eps << endl;
    if(eps < loweps){
      loweps = eps;
      lowest1 = tempvec1;
      lowest2 = tempvec2;
      lowest3 = tempvec3;
      lowalpha = alpha;
    }
  }//rotation
  //outfile2 << endl;
  
  cout << "Error:  " << setw(3) << low[0] << setw(12) << low[1] << setw(12) << low[2] << setw(12) << loweps << endl;
  cout << "Length: " << setw(3) << t1 << setw(12) << t2 << setw(12) << t3 << endl;
  cout << "Volume: " << setw(3) << abs(DotProduct(t1*n1, CrossProduct(t3*n3, t2*n2))) << endl;
  cout << endl;
  cout << setw(12) << "x" << setw(12) << "y" << setw(12) << "z" << setw(12) << "magnitude" << endl;
  cout << setw(12) << t1*n1.x << setw(12) << t1*n1.y << setw(12) << t1*n1.z << setw(12) << (t1*n1).Magnitude() << endl;
  cout << setw(12) << t2*n2.x << setw(12) << t2*n2.y << setw(12) << t2*n2.z << setw(12) << (t2*n2).Magnitude() << endl;
  cout << setw(12) << t3*n3.x << setw(12) << t3*n3.y << setw(12) << t3*n3.z << setw(12) << (t3*n3).Magnitude() << endl << endl;
  
  cout << "After rotation " << lowalpha*180/M_PI << endl;
  cout << setw(12) << lowest1.x << setw(12) << lowest1.y << setw(12) << lowest1.z << setw(12) << lowest1.Magnitude() << endl;
  cout << setw(12) << lowest2.x << setw(12) << lowest2.y << setw(12) << lowest2.z << setw(12) << lowest2.Magnitude() << endl;
  cout << setw(12) << lowest3.x << setw(12) << lowest3.y << setw(12) << lowest3.z << setw(12) << lowest3.Magnitude() << endl;
  cout << endl;

  return 0;
}



/*
void lattice_params(double a, double b, double c, double alpha, double beta, double gamma){
  Vector vec1 = Vector(a, 0 ,0);
  Vector vec2 = Vector(b*cos(gamma), b*sin(gamma), 0);
  Vector vec3 = Vector(
		       c*cos(beta),
		       c*(cos(alpha) - cos(gamma)*cos(beta))/sin(gamma),
		       sqrt(c*c - (c*cos(beta)*c*cos(beta)) - (c*(cos(alpha)-cos(gamma)*cos(beta))/sin(gamma)*c*(cos(alpha) - cos(gamma)*cos(beta))/sin(gamma)))
		       );

  
  
};
*/
