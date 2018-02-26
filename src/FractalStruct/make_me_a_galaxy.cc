#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void make_me_a_galaxy(int FractalRank,int numbers,double total_mass,double G,
			vector <double>& posx,vector <double>& posy,vector <double>& posz,
			vector <double>& velx,vector <double>& vely,vector <double>& velz)
  {
    int seed=9973+256*FractalRank;
    std::mt19937 genRAN(seed);
    std::uniform_real_distribution<double> distRAD(1.0e-10,1.0);
    std::uniform_real_distribution<double> distTHETA(-0.99999999,0.99999999);
    std::uniform_real_distribution<double> distPHI(0.0,8.0*atan(1.0));
    std::normal_distribution<double> distNORM(0.0,1.0);
    double rmax=30.0;
    double x_off=-1.0;
    double y_off=1.5;
    double z_off=0.5;
    double slope=-0.9;
    double velratio=0.7;
    double sigratio=0.002;
    velratio*=sqrt(G);
    sigratio*=sqrt(G);
    double slope3=slope+3.0;
    double expo=1.0/(3.0+slope);
    for(int ni=0;ni<numbers;ni++)
      {
	double r1=distRAD(genRAN);
	r1=pow(r1,expo);
	double ctheta=distTHETA(genRAN);
	double phi=distPHI(genRAN);
	double r=rmax*r1;
	posz[ni]=r*ctheta+z_off;
	double stheta=sqrt(1.0-ctheta*ctheta);
	double cphi=cos(phi);
	double sphi=sin(phi);
	posx[ni]=r*stheta*cphi+x_off;
	posy[ni]=r*stheta*sphi+y_off;
	double massr=total_mass*pow((r/rmax),slope3);
	double sqmr=sqrt(massr/r);
	double vt=velratio*sqmr;
	double vs=sigratio*sqmr;
	velx[ni]=-sphi*vt+vs*distNORM(genRAN);
	vely[ni]=cphi*vt+vs*distNORM(genRAN);
	velz[ni]=vs*distNORM(genRAN);
      }
  }
}
