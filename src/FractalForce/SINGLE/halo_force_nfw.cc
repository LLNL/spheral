#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  string Fractal::force_fixed="nfw";
  //
  double NFWmass(const double& r);
  double NFWforce(const double& r);
  double NFWpot(const double& r);
  void halo_force_fixed(Fractal& frac)
  {
    double NFW_rs=frac.get_halo_scale();
    double NFW_rho0=frac.get_halo_density0();
    double fourpi=16.0*atan(1.0);
    double NFW_pot_const=fourpi*NFW_rho0*NFW_rs*NFW_rs;
    double NFW_force_const=fourpi*NFW_rho0*NFW_rs;
    for (int n=0;n<frac.get_number_particles();n++)
      {
	vector <double> pos(3);
	vector <double> field(4);
	Particle* p=frac.particle_list[n];
	p->get_pos(pos);
	const double r=p->get_r(0.5,0.5,0.5);
	const double x=r/NFW_rs;
	double forcer=NFW_force_const*NFWforce(x)/r;
	field[0]=NFW_pot_const*NFWpot(x);
	field[1]=(pos[0]-0.5)*forcer;
	field[2]=(pos[1]-0.5)*forcer;
	field[3]=(pos[2]-0.5)*forcer;
	p->add_field_pf(field);
      }
  }
  inline double NFWmass(const double& x)
  {
    return log(1.0+x)-x/(1.0+x);
  }
  inline double NFWforce(const double& x)
  {
    return 1.0/((1.0+x)*x)-log(1.0+x)/(x*x);
  }
  inline double NFWpot(const double& x)
  {
    return -log(1.0+x)/x;
  }
}
