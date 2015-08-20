#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void make_me_some_particles(int rank,int numbers,double total_mass,vector <double>& masses,double G,
		     vector <double>& posx,vector <double>& posy,vector <double>& posz,
		     vector <double>& velx,vector <double>& vely,vector <double>& velz)
  {
    if(rank == 0)
      {
	posx[0]=0.2;
	posy[0]=0.3;
	posz[0]=0.4;
	velx[0]=0.0;
	velx[1]=0.0;
	velx[2]=0.0;
      }
    else
      {
	posx[0]=0.2;
	posy[0]=0.3;
	posz[0]=0.7;
	velx[0]=0.0;
	velx[1]=0.0;
	velx[2]=0.0;
      }
    /*
    double rand_max=(double)RAND_MAX;
    double rmax=10.0;
    double x_off=10.0;
    double y_off=-28.0;
    double z_off=11.0;
    double slope=-0.5;
    double velratio=0.5;
    velratio*=sqrt(G);
    double slope3=slope+3.0;
    double expo=1.0/(3.0+slope);
    double twopi=8.0*atan(1.0);
    for(int ni=0;ni<numbers;ni++)
      {
	double r1=Fractal::my_rand(rand_max);
	r1=pow(r1,expo);
	double ctheta=2.0*Fractal::my_rand(rand_max)-1.0;
	double phi=twopi*Fractal::my_rand(rand_max);
	double r=rmax*r1;
	posz[ni]=r*ctheta+z_off;
	double stheta=sqrt(1.0-ctheta*ctheta);
	double cphi=cos(phi);
	double sphi=sin(phi);
	posx[ni]=r*stheta*cphi+x_off;
	posy[ni]=r*stheta*sphi+y_off;
	double massr=total_mass*pow((r/rmax),slope3);
	double vt=velratio*sqrt(massr/r);
	velx[ni]=-sphi*vt;
	vely[ni]=cphi*vt;
	velz[ni]=0.0;
      }
    */
  }
}
