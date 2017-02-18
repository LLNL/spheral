#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
#include "interface_headers.hh"

//! Generic interface to FractalForce, used to be FractalGravity, Poisson solver.
//! The periodic version is not yet working. Patience.
//! 
//! The Fractal_Memory object contains all the necessary information.
//! It should be constructed before the first call to fractal_force and only deleted after last use.
//! The object is constructed with "fractal_interface_setup".
//! User must supply the 13 parameters needed.
//! 
//! The Fractal object is constructed with "fractal_create".
//! No user input needed.
//! 
//! It must be filled with Particles with "fractal_populate".
//! The user must supply the code to transfer positions and velocities to fractal.
//! 
//! Now you can run "fractal_force".
//! 
//! Then you call "fractal_depopulate" to transfer potential and forces back to your data structures.
//! The user must supply the code to transfer potential and forces to your data structures.
//! 
//! Then delete the Fractal object with "fractal_delete".
//! Do not delete the Fractal_Memory object.
//! 
//! For each step you do
//! 
//! fractal_create
//! fractal_populate
//! fractal_force
//! fractal_depopulate
//! fractal_delete
//! 
//! When you are finished call
//! fractal_memory_content_delete
//! 
//! End with 
//! delete Fractal_Memory object.
//! 
//! Parameters for fractal_interface_setup
//! periodic        ( true if periodic BC, false if isolate BC)
//! grid_length     ( length of fundamental grid, must be even, maximum 1024)
//! level_max       ( number of levels of refinement)
//! minimum_number  ( minimum number of particles belonging to a high point)
//! padding         (-1) => resolution jumps by factors of 2 (cheap)
//!                 (1)  => each high point fully padded (expensive)
//!                 (0)  => no padding (not recommended)
//! 
//! tolHypre        ( accuracy of hypre solver (1.0e-7 recommended))
//! maxitsHypre     ( max number of iterations in Hypre (20 recommended))
//! fractalDebug    ( debug run (true) extra IO)
//! FractalNodes0   ( Number of processors in x-direction)
//! FractalNodes1   ( Number of processors in y-direction)
//! FractalNodes2   ( Number of processors in z-direction)
//! BaseDirectory   ( (string) Base directory for Fractal output, absolute path)
//! RunIdentifier   ( (string) Unique Identifier, no spaces)
//! 
//! During simulation the following parameters can be changed, BEFORE fractal_create.
//! level_max, minimum_number, padding, tolHypre, maxitsHypre, fractalDebug
//! They are all public members of FractalMemory object.
//! For example, PFM->tolHypre=????;
//! When the fractal object is constructed, all parameters are transferred to the Fractal object
//! 
//! The call to fractal_force is padded with timing calls.
//! 
//!     PF->timing(-2,0);
//!     PF->timing(-1,49);
//!     fractal_force(*PF,*PFM);
//!     PF->timing(1,49);
//!     PF->timing(0,0);
//!     PF->timing_lev(0,0);
//! 
//! 
namespace FractalSpace
{
  Fractal_Memory* fractal_interface_setup(bool periodic,
					  int grid_length,
					  int level_max,
					  int minimum_number,
					  int padding,
					  int tolHypre,
					  int maxitsHypre,
					  bool fractalDebug,
					  int FractalNodes0,
					  int FractalNodes1,
					  int FractalNodes2,
					  string BaseDirectory,
					  string RunIdentifier)
    {
      // Construct a Fractal_Memory object. 
      // This will be used throughout the simulation.
      // All objects and information in fractal can be accessed from pointers stored in *PFM.
      Fractal_Memory* PFM= new Fractal_Memory;
      
      // Set up the parameters for Fractal. 
      fractal_memory_interface_parameters(*PFM,
					  periodic,
					  grid_length,
					  level_max,
					  minimum_number,
					  padding,
					  tolHypre,
					  maxitsHypre,
					  fractalDebug,
					  FractalNodes0,
					  FractalNodes1,
					  FractalNodes2,
					  BaseDirectory,
					  RunIdentifier);

      // Construct a Mess object. 
      // All MPI stuff is done in Mess member functions. 
      // Keeps all Box information
      // This will be used throughout the simulation.
      Mess* p_mess=new Mess(PFM->MPIrun,PFM->grid_length,PFM->periodic,PFM->number_particles);
      PFM->p_mess=p_mess;

      // Construct a File object. 
      // All output is done in File member functions. 
      // This will be used throughout the simulation.
      File* p_file=new File(PFM->BaseDirectory,p_mess->FractalRank,PFM->RUN);
      PFM->p_file=p_file;

      // Calculate all simulation information needed. 
      // Includes Boxes, FFTW startup etc.
      // This will be used throughout the simulation.
      fractal_force_init(PFM);

      // Return to sender
      return PFM;
  }
}
namespace FractalSpace
{
  void fractal_memory_interface_parameters(Fractal_Memory& mem,
					   bool periodic,
					   int grid_length,
					   int level_max,
					   int minimum_number,
					   int padding,
					   int tolHypre,
					   int maxitsHypre,
					   bool fractalDebug,
					   int FractalNodes0,
					   int FractalNodes1,
					   int FractalNodes2,
					   string BaseDirectory,
					   string RunIdentifier)
  {
    if(periodic)
      {
	assert(0);
      }
    else
      {
	// use this part if using isolated Boundary Conditions
	//
	/***********************/
	//You must set these parameters
	mem.BaseDirectory=BaseDirectory;
	//base directory
	mem.RUN=RunIdentifier;
	//directory name desriptor
	mem.FractalNodes0=FractalNodes0;
	//number of nodes in x-direction
	mem.FractalNodes1=FractalNodes1;
	//number of nodes in y-direction
	mem.FractalNodes2=FractalNodes2;
	//number of nodes in z-direction
	mem.FractalNodes=mem.FractalNodes0*mem.FractalNodes1*mem.FractalNodes2;
	//total number of nodes
	mem.MPIrun=mem.FractalNodes > 1;
	// makes sense, if I have more than one node it is an MPIrun
	mem.periodic = false ;
	//true for periodic BC and false for isolated BC
	mem.grid_length = grid_length ;    
	// length of fundamental grid, must be even
	//mem.number_particles = 1 ; 
	//    mem.level_max = 0 ; 
	mem.level_max = level_max ; 
	mem.global_level_max=mem.level_max;
	// highest level, nothing magical about this
	// total resolution = grid_length*2**level_max
	mem.minimum_number = minimum_number ; 
	// min number of particles in a cell to make it a high density cell
	// minimum_number**(1/3) is the local resolution in units of the local
	// mean interparticle spacing
	mem.padding = padding ;
	mem.padding=min(mem.padding,1);
	// if (0) no padding
	// if (-1) high cells are selectively padded so that resolution never jumps more than factor of 2.
	// if (padding = 1) each high cell is padded
	// padding has to be (-1,0 or 1)in an MPI version
	mem.maxits = maxitsHypre ;
	// maximum number of iterations in Poisson Solver
	mem.epsilon_sor = tolHypre ;
	//convergence criterion in Poisson Solver
	mem.debug=FractalDebug;
	// Does extra testing and prints out a bunch of diagnostics
	//
	/***********************/
	//Leave these parameters alone
	mem.min_hypre_group_size=1;
	//    mem.min_hypre_group_size=-1;
	// minum group size to use hypre
	mem.new_points_gen=9;
	//Generate this many Points in each go
	mem.steps=0;
	mem.momentum_conserve=false;
	mem.amnesia=true; // (true) forget everything after you are done. (false) remember everything.
	mem.mind_wipe=false; // (true) delete everything and then come back without calculating anything.
	mem.fixed_potential=false; // (true) use the fixed potential.
	mem.calc_shear=false;// (true) if we calculate shear of force field
	mem.calc_density_particle=false;
	mem.do_vel=false;
	mem.start_up=false;
	mem.halo_fixed=false;
	//
	mem.splits=0;
	mem.masks=0;
	//
	mem.masks_init=0;
      }
  }
}
namespace FractalSpace
{
  // User must write her own code to input all Particle information to the Fractal object, this means postions (3D) and mass.
  // This is template to show what I need. Substitute your own code.
  template <class YDa,class YDb> void fractal_populate(Fractal_Memory* PFM, 
						       const int number_particles,
						       YDa& YourDataa, 
						       YDb& YourDatab)
  {
    vector <double> xmin(3);
    vector <double> xmax(3);
    vector <double> pos(3);
    double mass;
    Fractal* PF=PFM->p_fractal;
    Particle* Parts_in=new Particle[number_particles];
    PFM->p_mess->parts_interface=Parts_in;
    PFM->number_particles=number_particles;
    PF->set_number_particles(number_particles);
    //
    YourNameSpace::get_xmin_xmax(YourDataa,YourDatab,xmin,xmax);
    //
    double dinv=1.0/(xmax[0]-xmin[0]);
    double total_mass=0.0;
    for(int ni=0;ni<number_particles;ni++)
      {
	//
	YourNameSpace::get_pos_mass(ni,YourDataa,YourDatab,pos,mass);
	//
	total_mass+=mass;
	pos[0]=(pos[0]-xmin[0])*dinv;
	pos[1]=(pos[1]-xmin[1])*dinv;
	pos[2]=(pos[2]-xmin[2])*dinv;
	Particle* P=&Parts_in[ni];
	PF->particle_list.push_back(P);
	P->set_pos(pos);
	P->set_mass(mass);
      }
    PFM->total_mass=total_mass;
  }
}
namespace FractalSpace
{
  template void fractal_populate(Fractal_Memory* PFM, 
				 const int number_particles,
				 Fractal& YourDataa, 
				 Fractal& YourDatab);
}
namespace FractalSpace
{
  // User must write her own code to input all Particle information to your data object(s), this means forces (3D) and potential.
  // This is template to show what I need. Substitute your own code.
  template <class YDa,class YDb> void fractal_depopulate(Fractal_Memory* PFM, 
							 YDa& YourDataa, YDb& YourDatab)
  {
    vector <double> xmin(3);
    vector <double> xmax(3);
    vector <double> potforce(4);
    Fractal* PF=PFM->p_fractal;
    int number_particles=PFM->number_particles;
    //
    YourNameSpace::get_xmin_xmax(YourDataa,YourDatab,xmin,xmax);
    //
    double dinv=1.0/(xmax[0]-xmin[0]);
    //
    double G=YourNameSpace::get_grav_constant();
    //
    double conv_pot=dinv*G;
    double conv_force=dinv*dinv*G;
    for(int ni=0;ni<number_particles;ni++)
      {
	Particle* P=PF->particle_list[ni];
	P->get_field_pf(potforce);
	potforce[0]*=conv_pot;
	potforce[1]*=conv_force;
	potforce[2]*=conv_force;
	potforce[3]*=conv_force;
	//
	YourNameSpace::set_pot_force(ni,YourDataa,YourDatab,potforce);
	//
      }
  }
}
namespace FractalSpace
{
  template void fractal_depopulate(Fractal_Memory* PFM, 
				   Fractal& YourDataa, 
				   Fractal& YourDatab);
}
namespace FractalSpace
{
  void fractal_delete(Fractal_Memory* PFM)
  {
    Fractal* PF=PFM->p_fractal;
    Particle* P=PFM->p_mess->parts_interface;
    delete [] P;
    P=0;
    delete PF;
    PF=0;
  }
}
namespace FractalSpace
{
  void fractal_create(Fractal_Memory* PFM)
  {
    Fractal* PF=new Fractal(*PFM);
    PFM->p_fractal=PF;
  }
}
namespace FractalSpace
{
  void fractal_memory_content_delete(Fractal_Memory* PFM)
  {
    Fractal* PF=PFM->p_fractal;
    Particle* P=PFM->p_mess->parts_interface;
    delete [] P;
    P=0;
    delete PF;
    PF=0;
    Mess* p_mess=PFM->p_mess;
    delete p_mess;
    p_mess=0;
    File* p_file=PFM->p_file;
    delete p_file;
    p_file=0;
  }
}
