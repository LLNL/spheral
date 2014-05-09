#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void make_sphere(Fractal& fractal)
  {
    double rand_max=(double)RAND_MAX;
    //    double ratio=fractal.get_length_ratio();
    srand(1234);
    double rmax=0.4;
    double x_off=0.51;
    double y_off=0.49;
    double z_off=0.52;
    double slope=-2.0;
    double mtot=0.5;
    //
    double expo=1.0/(3.0+slope);
    double m=mtot/fractal.get_number_particles();
    //  m=0.0;
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
	fractal.set_pos_x(n,x);
	fractal.set_pos_y(n,y);
	fractal.set_pos_z(n,z);
	fractal.set_particle_mass(n,m);
      }
    //   int n=0;
    //   fractal.set_pos_x(n,x_off);
    //   fractal.set_pos_y(n,y_off);
    //   fractal.set_pos_z(n,z_off); 
    //   fractal.set_particle_mass(n,mtot);
    //    cout << "m= " << m << "\n";
  }
}
