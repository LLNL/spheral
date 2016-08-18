#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void make_particles(Fractal_Memory& mem,Fractal& frac,int& count,const double& m,const bool& haha)
  {
    ofstream& FileFractal=fractal.p_file->FileFractal;
    double rand_max=(double)RAND_MAX;
    srand(1234);
    double r0=0.01;
    double rmax=0.4;
    double x_off=0.51;
    double y_off=0.49;
    double z_off=0.52;
    double sigma_e=0.0;
    int number=frac.get_number_particles();
    double mass=number*m;
    //
    double dens=-1.0;
    double x1=-1.0;
    double xmax=rmax/r0;
    double twopi=8.0*atan(1.0);
    //
    for(int n=0;n<number;++n)
      {
	bool test=false;
	while(!test)
	  {
	    x1=xmax*Frac::my_rand(rand_max);
	    dens=4.0*x1/(1.0+x1)/(1.0+x1);
	    test=Frac::my_rand(rand_max) <= dens;
	  }
	double ctheta=2.0*Frac::my_rand(rand_max)-1.0;
	double phi=twopi*Frac::my_rand(rand_max);
	double r=r0*r1;
	double z=r*ctheta+z_off;
	double stheta=sqrt(1.0-ctheta*ctheta);
	double x=r*stheta*cos(phi)+x_off;
	double y=r*stheta*sin(phi)+y_off;
	frac.particle_list[n]->set_pos(x,y,z);
	frac.particle_list[n]->set_mass(m);
      }
    FileFractal << " particle m= " << m << "\n";
  }
}
