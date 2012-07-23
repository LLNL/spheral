#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  string Fractal::sim_parameters="cosmos_flat_lambda";
  //
  template <class M> void fractal_memory_parameters(M& mem)
  {
    //    mem.padding=1; 
    mem.periodic=false;
    mem.debug=false;
    //    mem.number_steps_total=5;
    mem.new_points_gen=9;
    mem.remember_points=true;
    mem.number_steps_total=903;
    mem.redshift_start=99.0;
    mem.max_particles=300000;
    mem.omega_0=0.3;
    mem.omega_lambda=0.7;
    mem.h=0.7;
    mem.steps=-1;
    mem.random_gen=54321;
    // sets your values for parameters
    mem.amnesia=true; // (true) forget everything after you are done. (false) remember everything.
    mem.mind_wipe=false; // (true) delete everything and then come back without calculating anything.
    mem.fixed_potential=false; // (true) use the fixed potential.
    mem.calc_shear=false;// (true) if we calculate shear of force field
    mem.calc_density_particle=true;
    mem.do_vel=true;
    mem.start_up=true;
    mem.halo=mem.halo && !mem.periodic;
    mem.halo_fixed=mem.halo_fixed && !mem.periodic;
    mem.length_ratio=1;
    mem.omega_start=Omega(mem.omega_0,mem.omega_lambda,mem.redshift_start);
    mem.lambda_start=Lambda(mem.omega_0,mem.omega_lambda,mem.redshift_start);
    mem.sigma_initial=mem.sigma_0*Growth(mem.omega_0,mem.omega_lambda,mem.redshift_start);
    mem.time=Age_of_the_universe(mem.omega_start,mem.lambda_start,0.0);
    mem.total_mass=1.0;
    cout << " cosmo " << mem.omega_start << " " << mem.lambda_start << " " << mem.sigma_initial << " " << mem.time << endl;
    //
    mem.crash_levels=5;
    mem.crash_pow=2.0;
    mem.density_crash=5.5;
    mem.splits=2;
    mem.splits_center_x.assign(mem.splits,0.5);
    mem.splits_center_y.assign(mem.splits,0.5);
    mem.splits_center_z.assign(mem.splits,0.5);
    mem.splits_rad_x.assign(mem.splits,0.5);
    mem.splits_rad_y.assign(mem.splits,0.5);
    mem.splits_rad_z.assign(mem.splits,0.5);
    mem.splits_many.assign(mem.splits,1);
    mem.splits_square.assign(mem.splits,true);
    mem.splits_center_x[0]=0.05;
    mem.splits_center_y[0]=0.05;
    mem.splits_center_z[0]=0.05;
    mem.splits_rad_x[0]=0.1;
    mem.splits_rad_y[0]=0.1;
    mem.splits_rad_z[0]=0.1;
    mem.splits_many[0]=2;
    mem.splits_square[0]=false;
    mem.splits_center_x[1]=0.05;
    mem.splits_center_y[1]=0.05;
    mem.splits_center_z[1]=0.05;
    mem.splits_rad_x[1]=0.07;
    mem.splits_rad_y[1]=0.07;
    mem.splits_rad_z[1]=0.07;
    mem.splits_many[1]=3;
    mem.splits_square[1]=false;
    //
    mem.masks=4;
    mem.masks_center_x.assign(mem.masks,0.5);
    mem.masks_center_y.assign(mem.masks,0.5);
    mem.masks_center_z.assign(mem.masks,0.5);
    mem.masks_rad_x.assign(mem.masks,0.5);
    mem.masks_rad_y.assign(mem.masks,0.5);
    mem.masks_rad_z.assign(mem.masks,0.5);
    mem.masks_level.assign(mem.masks,0);
    mem.masks_square.assign(mem.masks,true);
    mem.masks_center_x[0]=0.5;
    mem.masks_center_y[0]=0.5;
    mem.masks_center_z[0]=0.5;
    mem.masks_rad_x[0]=0.5;
    mem.masks_rad_y[0]=0.5;
    mem.masks_rad_z[0]=0.5;
    mem.masks_level[0]=0;
    mem.masks_square[0]=true;
    mem.masks_center_x[1]=0.5;
    mem.masks_center_y[1]=0.5;
    mem.masks_center_z[1]=0.5;
    mem.masks_rad_x[1]=0.35;
    mem.masks_rad_y[1]=0.35;
    mem.masks_rad_z[1]=0.35;
    mem.masks_level[1]=2;
    mem.masks_square[2]=true;
    mem.masks_center_x[2]=0.5;
    mem.masks_center_y[2]=0.5;
    mem.masks_center_z[2]=0.5;
    mem.masks_rad_x[2]=0.25;
    mem.masks_rad_y[2]=0.25;
    mem.masks_rad_z[2]=0.25;
    mem.masks_level[2]=4;
    mem.masks_square[2]=true;
    mem.masks_center_x[3]=0.5;
    mem.masks_center_y[3]=0.5;
    mem.masks_center_z[3]=0.5;
    mem.masks_rad_x[3]=0.15;
    mem.masks_rad_y[3]=0.15;
    mem.masks_rad_z[3]=0.15;
    mem.masks_level[3]=8;
    mem.masks_square[3]=true;
  //
    mem.splits=0;
    //    mem.masks=0;
  //
    mem.masks_init=3;
    mem.masks_center_x_init.assign(mem.masks,0.5);
    mem.masks_center_y_init.assign(mem.masks,0.5);
    mem.masks_center_z_init.assign(mem.masks,0.5);
    mem.masks_rad_x_init.assign(mem.masks,0.5);
    mem.masks_rad_y_init.assign(mem.masks,0.5);
    mem.masks_rad_z_init.assign(mem.masks,0.5);
    mem.masks_level_init.assign(mem.masks,0);
    mem.masks_square_init.assign(mem.masks,true);
    mem.masks_center_x_init[0]=0.5;
    mem.masks_center_y_init[0]=0.5;
    mem.masks_center_z_init[0]=0.5;
    mem.masks_rad_x_init[0]=0.5;
    mem.masks_rad_y_init[0]=0.5;
    mem.masks_rad_z_init[0]=0.5;
    mem.masks_level_init[0]=0;
    mem.masks_square_init[0]=true;
    mem.masks_center_x_init[1]=0.5;
    mem.masks_center_y_init[1]=0.5;
    mem.masks_center_z[1]=0.5;
    mem.masks_rad_x_init[1]=0.3;
    mem.masks_rad_y_init[1]=0.3;
    mem.masks_rad_z_init[1]=0.3;
    mem.masks_level_init[1]=2;
    mem.masks_square_init[2]=true;
    mem.masks_center_x_init[2]=0.5;
    mem.masks_center_y_init[2]=0.5;
    mem.masks_center_z_init[2]=0.5;
    mem.masks_rad_x_init[2]=0.2;
    mem.masks_rad_y_init[2]=0.2;
    mem.masks_rad_z_init[2]=0.2;
    mem.masks_level_init[2]=4;
    mem.masks_square_init[2]=true;
  }
}
namespace FractalSpace
{
  template void  fractal_memory_parameters(Fractal_Memory& mem);
}
