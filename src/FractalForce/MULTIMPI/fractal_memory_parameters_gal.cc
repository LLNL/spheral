#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  string Fractal::sim_parameters="galaxy run";
  //
  template <class M> void fractal_memory_parameters(M& mem)
  {
    // These are the parameters you need to set.
    // For the others you should use my values for now.
    //    ofstream& FileFractal=mem.p_fractal->p_file->FileFractal;
    cout << " enter parameters " << endl;
    mem.BaseDirectory="/p/lscratchc/jensv/";
    //base directory
    mem.RUN="gal_cab_mpi";
    //directory name desriptor
    mem.MPIrun=true;
    // Is this an MPI run.
    mem.FractalNodes0=4;
    //number of nodes in x-direction
    mem.FractalNodes1=3;
    //number of nodes in y-direction
    mem.FractalNodes2=2;
    //number of nodes in z-direction
    mem.FractalNodes=mem.FractalNodes0*mem.FractalNodes1*mem.FractalNodes2;
    //total number of nodes
    cout << "FractalNodes " << mem.FractalNodes << endl;
    mem.MPIrun=mem.MPIrun || mem.FractalNodes > 1;
    // makes sense, if I have more than one node it is an MPIrun
    mem.periodic = false ;
    //true for periodic BC and false for isolated BC
    mem.max_particles=1000000;
    //The max number of particles the initial conditions code can generate through particle splitting
    //    mem.grid_length = 64 ;    
    mem.grid_length = 128 ;    
    // length of fundamental grid, must be even
    mem.number_particles = 90000 ; 
    //mem.number_particles = 1 ; 
    //    mem.level_max = 0 ; 
    mem.level_max = 8 ; 
    mem.global_level_max=mem.level_max;
    // highest level, nothing magical about this
    // total resolution = grid_length*2**level_max
    mem.minimum_number = 8 ; 
    // min number of particles in a cell to make it a high density cell
    // minimum_number**(1/3) is the local resolution in units of the local
    // mean interparticle spacing
    mem.padding = -1 ;
    mem.padding=min(mem.padding,1);
    // if (0) no padding
    // if (-1) high cells are selectively padded so that resolution never jumps more than factor of 2.
    // if (padding > 0) each high cell is padded by (2*padding+1)**3 cells.
    // padding has to be (-1,0 or 1)in an MPI version
    mem.min_hypre_group_size=5000;
    //    mem.min_hypre_group_size=-1;
    // minum group size to use hypre
    mem.maxits = 20 ;
    // maximum number of iterations in SOR algorithm, see NR
    mem.epsilon_sor = 1.0e-7 ;
    //convergence criterion in SOR, see NR
    mem.debug=true;
    // Does extra testing and prints out a bunch of diagnostics
    mem.new_points_gen=9;
    //Generate this many Points in each go
    mem.number_steps_total=400;
    // Total number of steps
    mem.number_steps_out=20;
    // Output how often
    mem.random_gen=54321;
    // Initial seed for std:rand
    // sets your values for parameters
    mem.steps=-1;
    mem.momentum_conserve=false;
    mem.amnesia=true; // (true) forget everything after you are done. (false) remember everything.
    mem.mind_wipe=false; // (true) delete everything and then come back without calculating anything.
    mem.fixed_potential=false; // (true) use the fixed potential.
    mem.calc_shear=false;// (true) if we calculate shear of force field
    mem.calc_density_particle=false;
    mem.do_vel=true;
    mem.start_up=true;
    mem.halo_fixed=mem.halo_fixed && !mem.periodic;
    mem.time=0.0;
    mem.step_length=0.0015;
    mem.total_mass=1.0;
    cout << " gal parameters " << endl;
    //
    mem.splits=0;
    mem.masks=0;
  //
    mem.masks_init=0;
  }
}
namespace FractalSpace
{
  template void fractal_memory_parameters(Fractal_Memory& mem);
}
