#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void make_particles(Fractal_Memory& fractal_memory,Fractal& fractal,int& count,const double& m,const bool& haha)
  {
    double rand_max=(double)RAND_MAX;
    srand(1234);
    double rmax=0.4;
    double x_off=0.51;
    double y_off=0.49;
    double z_off=0.52;
    double slope=-2.0;
    //
    double expo=1.0/(3.0+slope);
    double twopi=8.0*atan(1.0);
    //
    for(int n=0;n<fractal.get_number_particles();++n)
      {
	double r1=Fractal::my_rand(rand_max);
	r1=pow(r1,expo);
	double ctheta=2.0*Fractal::my_rand(rand_max)-1.0;
	double phi=twopi*Fractal::my_rand(rand_max);
	double r=rmax*r1;
	double z=r*ctheta+z_off;
	double stheta=sqrt(1.0-ctheta*ctheta);
	double x=r*stheta*cos(phi)+x_off;
	double y=r*stheta*sin(phi)+y_off;
	fractal.particle_list[n]->set_pos(x,y,z);
	fractal.particle_list[n]->set_mass(m);
      }
    cout << " particle m= " << m << endl;
  }
}
