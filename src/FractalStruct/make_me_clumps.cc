#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void make_me_a_galaxy(int FractalRank,int numbers,double total_mass,vector <double>& masses,double G,
		     vector <double>& posx,vector <double>& posy,vector <double>& posz,
		     vector <double>& velx,vector <double>& vely,vector <double>& velz)
  {
    int seed=9973+256*FractalRank;
    srand(seed);
    double rand_max=(double)RAND_MAX;
    double rmax=10.0;
    double x_off=1.0;
    double y_off=-2.0;
    double z_off=3.0;
    double slope=-0.65;
    double velratio=0.3333;
    velratio*=sqrt(G);
    double slope3=slope+3.0;
    double expo=1.0/(3.0+slope);
    double twopi=8.0*atan(1.0);
    int counter=0;
    double ax[8]={-10.0,9.0,-8.0,7.0,-9.0,10.0,-7.0,8.0};
    double ay[8]={-10.0,-10.0,9.0,9.0,-9.0,-9.0,10.0,10.0};
    double az[8]={-10.0,-9.0,-9.0,-10.0,9.0,10.0,10.0,9.0};
    for(int ni=0;ni<numbers;ni++)
      {
	double r1=Fractal::my_rand(rand_max);
	r1=pow(r1,expo);
	double ctheta=2.0*Fractal::my_rand(rand_max)-1.0;
	double phi=twopi*Fractal::my_rand(rand_max);
	double r=rmax*r1;
	posz[ni]=r*ctheta+z_off+az[counter];
	double stheta=sqrt(1.0-ctheta*ctheta);
	double cphi=cos(phi);
	double sphi=sin(phi);
	posx[ni]=r*stheta*cphi+x_off+ay[counter];
	posy[ni]=r*stheta*sphi+y_off+ax[counter];
	double massr=total_mass*pow((r/rmax),slope3);
	double vt=velratio*sqrt(massr/r);
	velx[ni]=-sphi*vt;
	vely[ni]=cphi*vt;
	velz[ni]=0.0;
	counter=(counter+1)%8;
      }
  }
}
