#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void make_me_a_disk(int FractalRank,double G,
		      int numberH,double massH,double scaleH,double flatH,double maxrH,
		      int numberD,double massD,double scalerD,double scalezD,double maxrD,
		      vector <double>& posx,vector <double>& posy,vector <double>& posz,
		      vector <double>& velx,vector <double>& vely,vector <double>& velz)
  {
    int seed=9973+256*FractalRank;
    srand(seed);
    //    std::default_random_engine generator(seed);
    //    std::uniform_real_distribution<double> distribution(0.0,1.0);
    double rand_max=(double)RAND_MAX;
    double x_off=0.0;
    double y_off=0.0;
    double z_off=0.0;
    double twopi=8.0*atan(1.0);
    bool allok=true;
    int ni=0;
    for(int nH=0;nH<numberH;nH++)
      {
	bool notOK=true;
	while(notOK)
	  {
	    double r=maxrH*Fractal::my_rand_not_zero(rand_max);
	    double dd=1.0/(1.0+scaleH*scaleH/(r*r));
	    notOK=Fractal::my_rand(rand_max) > dd;
	  }
	double ctheta=0.99999*(2.0*Fractal::my_rand(rand_max)-1.0);
	double phi=twopi*Fractal::my_rand(rand_max);
	posz[ni]=r*flatH*ctheta+z_off;
	double stheta=sqrt(1.0-ctheta*ctheta);
	double cphi=cos(phi);
	double sphi=sin(phi);
	posx[ni]=r*stheta*cphi+x_off;
	posy[ni]=r*stheta*sphi+y_off;
	velx[ni]=0.0;
	vely[ni]=0.0;
	velz[ni]=0.0;
	ni++;
      }
    double ddmax=exp(-1.0);
    for(int nD=0;nD<numberD;nD++)
      {
	bool notOK=true;
	while(notOK)
	  {
	    double r=maxrD*Fractal::my_rand_not_zero(rand_max);
	    double dd=r*exp(-d/scalerD);
	    notOK=ddmax*Fractal::my_rand(rand_max) > dd;
	  }
	double phi=twopi*Fractal::my_rand(rand_max);
	double cphi=cos(phi);
	double sphi=sin(phi);
	posx[ni]=r*cphi+x_off;
	posy[ni]=r*sphi+y_off;
	double z=0.99999*(2.0*Fractal::my_rand(rand_max)-1.0);
	posz[ni]=0.5*scalezD*log((1.0+x)/(1.0-x));
	velx[ni]=0.0;
	vely[ni]=0.0;
	velz[ni]=0.0;
	ni++;
      }
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
	velx[ni]=-sphi*vt+v2*cos(ang);
	vely[ni]=cphi*vt+v2*sin(ang);
	velz[ni]=vs*sqrt(-2.0*log(Fractal::my_rand_not_zero(rand_max)))*cos(twopi*Fractal::my_rand(rand_max));
	if(isfinite(posx[ni]) && isfinite(posy[ni]) && isfinite(posz[ni]) && isfinite(velx[ni]) && isfinite(vely[ni]) && isfinite(velz[ni]))
	  continue;
	cout << "BAD GALAXY " << FractalRank << " " << ni << " " << posx[ni] << " " << posy[ni] << " " << posz[ni];
	cout << " " << velx[ni] << " " << vely[ni] << " " << velz[ni];
	cout << " " << r1 << " " << phi << " " << ctheta << " " << v2 << " " << ang << endl;
	allok=false;
	//	masses[ni]=m;
      }
    assert(allok);
  }
}
