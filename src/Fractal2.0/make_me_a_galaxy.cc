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
    //    int seed=9973;
    srand(seed);
    //    std::default_random_engine generator(seed);
    //    std::uniform_real_distribution<double> distribution(0.0,1.0);
    double rand_max=(double)RAND_MAX;
    double rmax=30.0;
    double x_off=1.0;
    double y_off=-2.0;
    double z_off=3.0;
    double slope=-0.65;
    double velratio=0.5;
    double sigratio=0.5;
    velratio*=sqrt(G);
    sigratio*=sqrt(G);
    double slope3=slope+3.0;
    double expo=1.0/(3.0+slope);
    double twopi=8.0*atan(1.0);
    bool allok=true;
    //    double m=total_mass/static_cast<double>(numbers);
    for(int ni=0;ni<numbers;ni++)
      {
	//	double r1=std::distribution(generator);
	double r1=Fractal::my_rand_not_zero(rand_max);
	r1=pow(r1,expo);
	//	double ctheta=2.0*std::distribution(generator)-1.0;
	double ctheta=0.99999*(2.0*Fractal::my_rand(rand_max)-1.0);
	//	double phi=twopi*std::distribution(generator);
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
	double vs=sigratio*sqrt(massr/r);
	double v2=vs*sqrt(-2.0*log(Fractal::my_rand_not_zero(rand_max)));
	double ang=twopi*Fractal::my_rand(rand_max);
	//	posx[ni]=-25.0+0.01;  //
	//	posy[ni]=0.0-0.02;   //
	//	posz[ni]=25.0+0.001;  //
	velx[ni]=-sphi*vt+v2*cos(ang);
	vely[ni]=cphi*vt+v2*sin(ang);
	velz[ni]=vs*sqrt(-2.0*log(Fractal::my_rand_not_zero(rand_max)))*cos(twopi*Fractal::my_rand(rand_max));
	if(isfinite(posx[ni]) && isfinite(posy[ni]) && isfinite(posz[ni]) && isfinite(velx[ni]) && isfinite(vely[ni]) && isfinite(velz[ni]))
	  continue;
	cerr << "BAD GALAXY " << FractalRank << " " << ni << " " << posx[ni] << " " << posy[ni] << " " << posz[ni];
	cerr << " " << velx[ni] << " " << vely[ni] << " " << velz[ni];
	cerr << " " << r1 << " " << phi << " " << ctheta << " " << v2 << " " << ang << endl;
	allok=false;
	//	masses[ni]=m;
      }
    assert(allok);
  }
}
