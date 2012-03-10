#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  string Fractal::force_fixed="flat";
  //
  double FLATmass(const double& r);
  double FLATforce(const double& r);
  double FLATpot(const double& r);
  void halo_force_fixed(Fractal& frac)
  {
    double FLAT_rs=frac.get_halo_scale();
    double FLAT_rho0=frac.get_halo_density0();
    double fourpi=16.0*atan(1.0);
    double FLAT_pot_const=fourpi*FLAT_rho0*FLAT_rs*FLAT_rs;
    double FLAT_force_const=fourpi*FLAT_rho0*FLAT_rs;
    for (int n=0;n<frac.get_number_particles();n++)
      {
	vector <double> pos(3);
	vector <double> field(4);
	Particle* p=frac.particle_list[n];
	p->get_pos(pos);
	const double r=p->get_r(0.5,0.5,0.5);
	const double x=r/FLAT_rs;
	double forcer=FLAT_force_const*FLATforce(x)/r;
	field[0]=FLAT_pot_const*FLATpot(x);
	field[1]=(pos[0]-0.5)*forcer;
	field[2]=(pos[1]-0.5)*forcer;
	field[3]=(pos[2]-0.5)*forcer;
	p->add_field_pf(field);
      }
  }
  inline double FLATmass(const double& x)
  {
    return x-atan(x);
  }
  inline double FLATforce(const double& x)
  {
    return -(x-atan(x))/x/x;
  }
  inline double FLATpot(const double& x)
  {
    return 0.5*log(1.0+x*x)+atan(x)/x-1.0;
  }
}
