#include "libs.hh"
#include "classes.hh"
#include "nbody.hh"
#include "headers.hh"
namespace FractalSpace
{
template <class M>  void energies(M& mem, Fractal& fractal);
template <class M>  void take_a_step(M& mem, Fractal& fractal);
  template <class M>  void split_particle(M& mem,const double& x0,const double& y0,const double& z0,int& count,const double& m,const int& n);
}
int main()
{
  using namespace FractalSpace;
  static ofstream FilePos;
  if(!FilePos.is_open())
    FilePos.open("pos.d");
  cout << "starting out " << endl;
  Nbody* p_nbody=new (nothrow) Nbody;
  assert(p_nbody);
  Nbody& nbody=*p_nbody;
  cout << "nbody created " << endl;
  nbody.omega_0=0.3;
  nbody.omega_lambda=0.7;
  nbody.h=0.7;
  nbody.steps=-1;

  // sets your values for parameters
  bool amnesia=true; // (true) forget everything after you are done. (false) remember everything.
  bool mind_wipe=false; // (true) delete everything and then come back without calculating anything.
  bool fixed_potential=false; // (true) use the fixed potential.
  bool calc_shear=false;// (true) if we calculate shear of force field
  int length_ratio=1;

  //  double hubble_start=Hubble(nbody.omega_0,nbody.omega_lambda,nbody.redshift_start);
  nbody.omega_start=Omega(nbody.omega_0,nbody.omega_lambda,nbody.redshift_start);
  nbody.lambda_start=Lambda(nbody.omega_0,nbody.omega_lambda,nbody.redshift_start);
  double sigma_start=nbody.sigma_0*Growth(nbody.omega_0,nbody.omega_lambda,nbody.redshift_start);
  
  nbody.time=Age_of_the_universe(nbody.omega_0,nbody.omega_lambda,nbody.redshift_start);

  //populate positions and masses
  double pi=atan(1.0)*4.0;
  double rand_max=(double)RAND_MAX;
  double m=0.375*nbody.omega_start/pi/(double)nbody.number_particles;
  //  nbody.resize_vectors(nbody.number_particles);
  int length=nbody.grid_length;
  int length_2=length*length;
  double delta=1.0/(double)length;
  srand(1234);
  int count=0;
  for(int n=0;n<nbody.number_particles;++n)
    {
      double x0,y0,z0;
      if(nbody.random_initial)
	{
	  x0=(double)(rand())/rand_max;
	  y0=(double)(rand())/rand_max;
	  z0=(double)(rand())/rand_max;
	}
      else
	{
	  x0=((double)(n % length)+nbody.off_x)*delta;
	  y0=((double)(n/length % length)+nbody.off_y)*delta;
	  z0=((double)(n/length_2)+nbody.off_z)*delta;
	}
      split_particle(nbody,x0,y0,z0,count,m,0);
    }
  cout << "size " << nbody.pos_x.size() << endl;
  nbody.number_particles=count;
  //  assert(0);
  //   for(int i=0;i<nbody.number_particles;i++)
  //     {
  //       cout << "ok a " << i << " " << nbody.pos_x[i] << " " << nbody.pos_y[i] << " " << nbody.pos_z[i] << endl;
  //     }
  //  nbody.resize_vectors(count);
  nbody.pos_x.resize(count);
  nbody.pos_y.resize(count);
  nbody.pos_z.resize(count);
  nbody.vel_x.resize(count,0.0);
  nbody.vel_y.resize(count,0.0);
  nbody.vel_z.resize(count,0.0);
  cout << "make parts " << count << endl;
  //   for(int i=0;i<nbody.number_particles;i++)
  //     {
  //       cout << "ok b " << i << " " << nbody.pos_x[i] << " " << nbody.pos_y[i] << " " << nbody.pos_z[i] << endl;
  //     }
  //call the gravity solver, with some timing calls
  Fractal* p_fractal=new Fractal;
  Fractal& fractal=*p_fractal;
  Fractal_Memory* p_fractal_memory= new Fractal_Memory;
  Fractal_Memory& fractal_memory=*p_fractal_memory;
  fractal_memory.set_generated_from(p_fractal);
  //
  fractal.set_number_particles(nbody.number_particles);
  fractal.set_grid_length(nbody.grid_length); 
  fractal.set_moat_0(nbody.moat_0);
  fractal.set_periodic(nbody.periodic);
  fractal.set_minimum_number(nbody.minimum_number);
  fractal.set_level_max(nbody.level_max);
  fractal.set_padding(nbody.padding);
  fractal.set_tweaks(nbody.tweaks);
  fractal.set_force_max(nbody.force_max);
  fractal.set_length_ratio(length_ratio);
  fractal_memory.set_amnesia(amnesia);
  fractal_memory.set_mind_wipe(mind_wipe);
  fractal_memory.set_fixed_potential(fixed_potential);
  fractal_memory.set_calc_shear(calc_shear);
  fractal.set_parameters_bool(512,true);      // first time power spectrum
  fractal.set_parameters_bool(513,true);      // startup
  fractal.set_parameters_bool(514,false);  // Gauss Filter
  fractal.set_parameters_bool(515,true);  // Do Variance
  fractal.set_parameters_bool(520,false); // do density at particle
  fractal.set_parameters_int(512,nbody.spectrum_number);    // spectrum number
  fractal.set_parameters_int(513,nbody.highest_level_init);    // highest_level_init
  fractal.set_parameters_int(514,nbody.norm_what);    // normalize what? rho(0), vel(1), force(2), force(3)
  fractal.set_parameters_double(512,nbody.power_slope); // power slope
  fractal.set_parameters_double(513,1.0e1); // power spectrum cutoff, not norm scale
    //  fractal.set_parameters_double(513,1.0e6) // power spectrum cutoff, not norm scale
  fractal.set_parameters_double(514,nbody.norm_scale); // norm scale
  fractal.set_parameters_double(515,sigma_start); // /sigma_\rho
  double delta_z=Growth(nbody.omega_0,nbody.omega_lambda,nbody.redshift_start);
  double vfratio=1.0/(1.5*nbody.omega_start)*dGrowthdT(nbody.omega_start,nbody.lambda_start,0.0);
  double omega_fraction=1.0/(1.5*nbody.omega_start);
  fractal.omega_fraction=omega_fraction;
  double scaling=1.0;
  if(nbody.spectrum_number==1)
    {
      double a1=pow(46.9*nbody.omega_0*nbody.h*nbody.h,0.67)*(1.0+pow(32.1*nbody.omega_0*nbody.h*nbody.h,-0.532));
      double a2=pow(12.0*nbody.omega_0*nbody.h*nbody.h,0.424)*(1.0+pow(45.0*nbody.omega_0*nbody.h*nbody.h,-0.582));
      double alpha=pow(a1,-nbody.omega_b/nbody.omega_0)*pow(a2,-pow(nbody.omega_b/nbody.omega_0,3));
      scaling=nbody.box_length*nbody.omega_0*nbody.h*nbody.h*sqrt(alpha)*
	pow(1.0-nbody.omega_b/nbody.omega_0,0.6);
      cout << "scaling " << a1 << " " << a2 << " " << alpha << " " << " " << nbody.box_length << " " << nbody.h << " " << scaling << endl;
    }
  //  assert(0);
  fractal.set_parameters_double(516,scaling);
  //set vector sizes, sets force and potential to zero for all particles
  fractal.resize_all();
  for(int n=0;n<nbody.number_particles;++n)
    {
      fractal.set_pos(n,nbody.pos_x[n],nbody.pos_y[n],nbody.pos_z[n]);
      fractal.set_particle_mass(n,nbody.particle_mass[n]);
    }
  double drad=(1.0+nbody.redshift_start)/(double)fractal.growth_length;
  for(int i=0;i<=fractal.growth_length;i++)
    {
      fractal.rad[i]=i*drad+1.0;
      fractal.grow[i]=Growth(nbody.omega_start,nbody.lambda_start,1.0/fractal.rad[i]-1.0);
      cout << i << " " << fractal.rad[i] << " " << fractal.grow[i] << endl;
    }
  //  assert(0);
  //
  fractal.timing(-2,0);
  fractal.timing(-1,29);
  fractal_gravity(fractal,fractal_memory);
  fractal.timing(1,29);
  fractal.timing(0,0);
  
  cout << "initial " << delta_z << " " << vfratio << " " << omega_fraction << endl;
  if(nbody.norm_what == 1)omega_fraction/=vfratio;
  for(int n=0;n<nbody.number_particles;++n)
    {
      fractal.add_pos_x(n,fractal.get_force_x(n)*omega_fraction);
      fractal.add_pos_y(n,fractal.get_force_y(n)*omega_fraction);
      fractal.add_pos_z(n,fractal.get_force_z(n)*omega_fraction);
    }
  fractal.set_parameters_bool(513,false);
  fractal.timing(-2,0);
  fractal.timing(-1,29);
  fractal_gravity(fractal,fractal_memory);
  fractal.timing(1,29);
  fractal.timing(0,0);
  double lambda=Lambda(nbody.omega_0,nbody.omega_lambda,nbody.redshift_start);
  double omega_fraction_v=dGrowthdT(nbody.omega_start,lambda,1.0/(1.0+nbody.step_length*0.5)-1.0);
  for(int n=0;n<nbody.number_particles;++n)
    {
      nbody.pos_x[n]+=fractal.get_force_x(n)*omega_fraction;
      nbody.pos_y[n]+=fractal.get_force_y(n)*omega_fraction;
      nbody.pos_z[n]+=fractal.get_force_z(n)*omega_fraction;
      nbody.vel_x[n]=fractal.get_force_x(n)*omega_fraction_v;
      nbody.vel_y[n]=fractal.get_force_y(n)*omega_fraction_v;
      nbody.vel_z[n]=fractal.get_force_z(n)*omega_fraction_v;
    }
  nbody.arad=1.0;
  nbody.time=Age_of_the_universe(nbody.omega_0,nbody.omega_lambda,nbody.redshift_start);
  fractal.set_parameters_bool(513,false);      // not startup
  //forget everything
  for(int step=0;step < nbody.number_steps_total; ++step)
    {
      fractal.set_parameters_bool(515,step % 10 == 9);  // Do Variance
      for(int n=0;n<nbody.number_particles;++n)
	{
	  fractal.set_pos_x(n,nbody.pos_x[n]);
	  fractal.set_pos_y(n,nbody.pos_y[n]);
	  fractal.set_pos_z(n,nbody.pos_z[n]);
	}
      if(nbody.steps == -1) energies(nbody,fractal);
      nbody.steps=step;
      cout << "step = " << step << " " << nbody.arad << endl;
      if(step % nbody.number_steps_out == 0)
	fractal.set_parameters_bool(520,true);
      else
	{
	  fractal.set_parameters_bool(520,false);
	  fractal.resize_dens(0);
	}
      fractal.timing(-2,0);
      fractal.timing(-1,29);
      fractal_gravity(fractal,fractal_memory);
      fractal.timing(1,29);
      fractal.timing(0,0);
      take_a_step(nbody,fractal);
      energies(nbody,fractal);
      if(step % nbody.number_steps_out == 0)
	{
	  for(int n=0;n<nbody.number_particles;++n)
	    {
	      FilePos << "stepout " << step << " " << n << " " << nbody.pos_x[n] << " " << nbody.pos_y[n] << " " << nbody.pos_z[n] ;
	      FilePos << " " << nbody.vel_x[n] << " " << nbody.vel_y[n] << " " << nbody.vel_z[n] << " " << nbody.particle_mass[n] ;
	      //		 FilePos << " lev= " << nbody.highest_level[n] << " ";
	      if(fractal.get_parameters_bool(520)) FilePos << " " << fractal.get_density(n);
	      FilePos  << endl;
	    }
	}
    }
  delete p_fractal;
  p_fractal=0;
  return 1;
}
