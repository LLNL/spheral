#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void make_sphere(Fractal& fractal);
}
int main()
{
  using namespace FractalSpace;
  /*
    This is a pseudo calling routine for the poisson solver:

    The calling interface is simple, you call the routine fractal_gravity
    with the Fractal class object "fractal" or whatever you want to call it.
    All real numbers are "double" and all integers are "int".

    void fractal_gravity(Fractal& fractal);

    "fractal" contains all the information the poisson solver needs and nothing else.
    The internal workings and data are hidden from the user as decreed by the creator.

    Input is positions and masses for all particles plus parameters.
    Output is potential and accelerations plus timings. Everything else is unchanged.
    "fractal" is created using the default constructor before each use and must be deleted after each use.

    Construct fractal
    Set the parameters below to your own values
    Do not touch other parameters in fractal
    Set vector sizes and null the vectors with fractal.resize_all()
    Assign the positions and masses
    Call fractal_gravity
    Read off forces(i.e. accelerations) and potential
    delete fractal

    For cosmo, at beginning of simulation (arad=1), Omega has value omega_0.
    The total mass in the simulations  is M_jvv=3*omega_0/(8*pi)
    The masses must be rescaled by the factor M_jvv/M_mike
    The gravitational constant is G_jvv=1
    The positions have to fit into the unit cube.
    x_jvv=(x_mike-x_mike_0)/scale_mike
    m_jvv=m_mike*(M_jvv/M_mike)

    So on output
    potential_mike=potential_jvv*G_mike*M_mike/M_jvv/scale_mike
    force_mike=force_jvv*G_mike*M_mike/M_jvv/scale_mike**2
  */
  int haha0=clock();
  Fractal* p_fractal=new Fractal;
  Fractal& fractal=*p_fractal;
  Fractal_Memory* p_fractal_memory= new Fractal_Memory;
  Fractal_Memory& fractal_memory=*p_fractal_memory;
  fractal_memory.set_generated_from(p_fractal);
  int haha1=clock();
  //
  // default values
  // replace with your own values
  int number_particles=262144;  // total numbe of particles
  int grid_length=64;   // length of grid for FFT solver, power of 2
  bool periodic=true;  // true for cosmological simulation
  //  bool periodic=true;  // true for cosmological simulation
  int minimum_number=8; // minimum number of particles in a highcell
  int level_max=8; // how many levels in the tree
  int padding=0; // (=0) no padding, (=1) full padding
  int tweaks=2; // (=2) making sure no interface overlaps. (=0) I don't care
  double epsilon_sor=6.0e-5; // accuracy parameter in SOR, see NR
  int maxits=250; // maximum number of iterations in SOR, see NR
  bool density_at_particle=false; // if true, calculate density at particle positions
  bool amnesia=false; // forget everything after you are done
  bool mind_wipe=false; // delete everythign and then come back
  bool fixed_potential=false; // use the fixed potential
  int length_ratio=4;
  //memory_value=0 on output if all groups converge.
  // if one of the groups does not convergence, then memory_value=8 on output.
  // sets your values for parameters
  fractal.set_number_particles(number_particles);
  fractal.set_grid_length(grid_length); 
  fractal.set_periodic(periodic);
  fractal.set_minimum_number(minimum_number);
  fractal.set_level_max(level_max);
  fractal.set_padding(padding);
  fractal.set_tweaks(tweaks);
  fractal.set_epsilon_sor(epsilon_sor);
  fractal.set_maxits(maxits);
  fractal.set_parameters_bool(520,density_at_particle);
  fractal.set_force_max(-10.0);
  fractal.set_length_ratio(length_ratio);
  fractal_memory.set_amnesia(amnesia);
  fractal_memory.set_mind_wipe(mind_wipe);
  fractal_memory.set_fixed_potential(fixed_potential);
  //set vector sizes, sets force and potential to zero for all particles
  fractal.resize_all();
  int haha3=clock();
  //populate positions and masses
  make_sphere(fractal);
  //call the gravity solver, with some timing calls
  fractal.timing(-2,0);
  fractal.timing(-1,29);
  fractal_gravity(fractal,fractal_memory);
  fractal.timing(1,29);
  fractal.timing(0,0);
  double ratio=length_ratio;
  for(int n=0;n<fractal.get_number_particles();++n)
    {
      double x=fractal.get_pos_x(n)-0.51;
      double y=fractal.get_pos_y(n)-0.49;
      double z=fractal.get_pos_z(n)-0.52/ratio;
      double r=sqrt(x*x+y*y+z*z);
      double pot_jvv=fractal.get_potential(n);
      double f_x=fractal.get_force_x(n);
      double f_y=fractal.get_force_y(n);
      double f_z=fractal.get_force_z(n);
      double fr=(f_x*x+f_y*y+f_z*z)/r;
      double fr_abs=abs(fr);
      double ftx=f_x-fr*x/r;
      double fty=f_y-fr*y/r;
      double ftz=f_z-fr*z/r;
      double ft=sqrt(ftx*ftx+fty*fty+ftz*ftz);
      if(n < 10000)
	{
	  cout << "resa " << n << " " << r << " " << pot_jvv <<" " << fr_abs <<  " " << ft << " " << fractal.get_particle_mass(n) << " ";
	  if(fractal.get_parameters_bool(520)) cout << fractal.get_density(n);
	  cout << " " << x <<  " " << y << " " << z ;
	  cout << endl;
	}
    }
  fractal_memory.set_fixed_potential(true);
  Fractal* p_fractal_asynch=new Fractal(fractal);
  Fractal& fractal_asynch=*p_fractal_asynch;
  fractal_asynch.set_number_particles(10000);
  fractal_asynch.resize_pospf();
  fractal_asynch.copy_pos(fractal,10000,20000);
  fractal_asynch.timing(-2,0);
  fractal_asynch.timing(-1,29);
  fractal_gravity(fractal_asynch,fractal_memory);
  fractal_asynch.timing(1,29);
  fractal_asynch.timing(0,0);
  for(int n=0;n<fractal_asynch.get_number_particles();++n)
    {
      double x=fractal_asynch.get_pos_x(n)-0.51;
      double y=fractal_asynch.get_pos_y(n)-0.49;
      double z=fractal_asynch.get_pos_z(n)-0.52/ratio;
      double r=sqrt(x*x+y*y+z*z);
      double pot_jvv=fractal_asynch.get_potential(n);
      double f_x=fractal_asynch.get_force_x(n);
      double f_y=fractal_asynch.get_force_y(n);
      double f_z=fractal_asynch.get_force_z(n);
      double fr=(f_x*x+f_y*y+f_z*z)/r;
      double fr_abs=abs(fr);
      double ftx=f_x-fr*x/r;
      double fty=f_y-fr*y/r;
      double ftz=f_z-fr*z/r;
      double ft=sqrt(ftx*ftx+fty*fty+ftz*ftz);
      cout << "resb " << n << " " << r << " " << pot_jvv <<" " << fr_abs <<  " " << ft << " " << endl;
    }
  int haha2=clock();
  delete p_fractal;
  p_fractal=0;
  delete p_fractal_asynch;
  p_fractal_asynch=0;
  int haha4=clock();
  cout << "fractal timing " << haha1-haha0 << " " << haha3-haha1 << " " << haha4-haha2 << endl;
  return 1;
}
