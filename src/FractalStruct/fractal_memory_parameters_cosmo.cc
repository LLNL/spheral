#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  string Fractal::sim_parameters="cosmos_flat_lambda";
  //
  void fractal_memory_parameters(Fractal_Memory* pmem,double _mulT_)
  {
    Fractal_Memory& mem=*pmem;
    int Ranky;
    MPI_Comm_rank(Fractal_Memory::FRACTAL_UNIVERSE,&Ranky);
    mem.FractalNodes=mem.FractalNodes0*mem.FractalNodes1*mem.FractalNodes2;
    // These are the parameters you need to set.
    // For the others you should use my values for now.
    //    ofstream& FileFractal=mem.p_fractal->p_file->FileFractal;
    //directory name descriptor
    // mem.MPIrun=true;
    mem.FFTNodes=9876543;
    mem.FFTNodes=min(mem.FFTNodes,mem.FractalNodes);
    mem.FFTNodes=min(mem.FFTNodes,mem.grid_length/2);
    // max number of nodes for FFTW
    mem.periodic = true ;
    //true for periodic BC and false for isolated BC
    if(mem.grid_length < 64)
      mem.grid_length = 256;
    // length of fundamental grid, must be even
    long int really_longGL=mem.grid_length;
    long int really_longFN=mem.FractalNodes;
    long int really_longGLFN=(really_longGL*really_longGL*really_longGL)/really_longFN;
    mem.number_particles=really_longGLFN;
    //    mem.number_particles = (mem.grid_length*mem.grid_length*mem.grid_length)/mem.FractalNodes; 
    // I will let you guess, you are wrong. It needs to be = grid_length**3/FractalNodes.
    mem.max_particles=(double)(mem.number_particles)*_mulT_+0.5;
    //The max number of particles the initial conditions code can generate through particle splitting
    mem.level_max = 10 ; 
    //    mem.level_max = 0 ; 
    mem.global_level_max=mem.level_max;
    mem.highest_level_init=min(mem.highest_level_init,mem.level_max);
    // highest level, nothing magical about this
    // total resolution = grid_length*2**level_max
    mem.minimum_number = 8 ; 
    // min number of particles in a cell to make it a high density cell
    // minimum_number**(1/3) is the local resolution in units of the local
    // mean interparticle spacing
    mem.padding = -1 ;
    // mem.padding = 1 ;
    mem.padding=min(mem.padding,1);
    // if (-1) high cells are selectively padded so that resolution never jumps more than factor of 2.
    // if (padding > 0) each high cell is padded by (2*padding+1)**3 cells.
    // padding has to be (0) or(1)in an MPI version
    mem.maxits = 20 ;
    // maximum number of iterations in SOR or Hypre Solver
    mem.HTOL = 1.0e-7;
    //convergence criterion in Hypre Solver
    mem.hypre_load_balance=true;
    //    mem.hypre_load_balance=false;
    mem.debug=true;
    // Does extra testing and prints out a bunch of diagnostics
    mem.new_points_gen=9;
    //Generate this many Points in each go
    mem.number_steps_total=906;
    //    mem.number_steps_total=8;
    // Total number of steps
    mem.number_steps_out=20;
    // Output how often
    mem.time_trial=true;
    mem.redshift_start=99.0;
    // step length
    mem.step_length=0.025;
    // initial redshift
    mem.omega_0=0.3;
    // Omega_matter at current epoch, NOT initial epoch.
    mem.omega_lambda=0.7;
    // Omega_Lambda at current epoch, NOT initial epoch.
    mem.h=0.7;
    // Hubble Constant in units of 100km/sec/Mpc at current epoch, NOT initial epoch.
    mem.random_gen=54321;
    // Initial seed for std:rand
    // sets your values for parameters
    mem.steps=-1;
    mem.min_hypre_group_size=729;
    mem.momentum_conserve=false;
    mem.amnesia=true; // (true) forget everything after you are done. (false) remember everything.
    mem.mind_wipe=false; // (true) delete everything and then come back without calculating anything.
    mem.fixed_potential=false; // (true) use the fixed potential.
    mem.calc_shear=true;// (true) if we calculate shear of force field
    mem.calc_density_particle=true;
    mem.do_vel=true;
    mem.start_up=true;
    mem.halo_fixed=mem.halo_fixed && !mem.periodic;
    mem.omega_start=Omega(mem.omega_0,mem.omega_lambda,mem.redshift_start);
    mem.lambda_start=Lambda(mem.omega_0,mem.omega_lambda,mem.redshift_start);
    mem.sigma_initial=mem.sigma_0*Growth(mem.omega_0,mem.omega_lambda,mem.redshift_start);
    mem.time=Age_of_the_universe(mem.omega_start,mem.lambda_start,0.0);
    mem.total_mass=1.0;
    if(Ranky == 0)
      cerr << " cosmo " << mem.omega_start << " " << mem.lambda_start << " " << mem.sigma_initial << " " << mem.time << "\n";
    //
    mem.crash_levels=8;
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
    mem.masks_rad_x[0]=1.5;
    mem.masks_rad_y[0]=1.5;
    mem.masks_rad_z[0]=1.5;
    mem.masks_level[0]=0;
    mem.masks_square[0]=true;
    mem.masks_center_x[1]=0.5;
    mem.masks_center_y[1]=0.5;
    mem.masks_center_z[1]=0.5;
    mem.masks_rad_x[1]=0.48;
    mem.masks_rad_y[1]=0.48;
    mem.masks_rad_z[1]=0.48;
    mem.masks_level[1]=4;
    mem.masks_square[2]=true;
    mem.masks_center_x[2]=0.5;
    mem.masks_center_y[2]=0.5;
    mem.masks_center_z[2]=0.5;
    mem.masks_rad_x[2]=0.46;
    mem.masks_rad_y[2]=0.46;
    mem.masks_rad_z[2]=0.46;
    mem.masks_level[2]=6;
    mem.masks_square[2]=true;
    mem.masks_center_x[3]=0.5;
    mem.masks_center_y[3]=0.5;
    mem.masks_center_z[3]=0.5;
    mem.masks_rad_x[3]=0.44;
    mem.masks_rad_y[3]=0.44;
    mem.masks_rad_z[3]=0.44;
    mem.masks_level[3]=10;
    mem.masks_square[3]=true;
  //
    mem.splits=0;
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
    mem.masks_center_z_init[1]=0.5;
    mem.masks_rad_x_init[1]=0.3;
    mem.masks_rad_y_init[1]=0.3;
    mem.masks_rad_z_init[1]=0.3;
    mem.masks_level_init[1]=2;
    mem.masks_square_init[1]=true;
    mem.masks_center_x_init[2]=0.5;
    mem.masks_center_y_init[2]=0.5;
    mem.masks_center_z_init[2]=0.5;
    mem.masks_rad_x_init[2]=0.2;
    mem.masks_rad_y_init[2]=0.2;
    mem.masks_rad_z_init[2]=0.2;
    mem.masks_level_init[2]=4;
    mem.masks_square_init[2]=true;
    mem.masks_init=0;

    if(Ranky == 0)
      cerr << " finishing cosmo " << "\n";
  }
}
