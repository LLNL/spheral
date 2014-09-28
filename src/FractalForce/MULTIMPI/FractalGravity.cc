#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
#include "fractal_interface_public.hh"
FractalSpace::Fractal_Memory* FractalGravityCreate();
    //! Parameters for FractalGravity interface_setup, MUST BE SET.
    //! JVV stuff
    //! NumberParticles ( Number of particles on this node only).
    //! balance         (0) => equal volume domains (easy but slow). Use for periodic BC for now.
    //!                 (1) => equal # particles + level 0 points for each domain (significant speedup for Isolated BC)
    //! FractalNodes0   ( Number of processors in x-direction)
    //! FractalNodes1   ( Number of processors in y-direction)
    //! FractalNodes2   ( Number of processors in z-direction)
    //!                 ( For Periodic we must have FractalNodes0=FractalNodes1=FractalNodes2; will be fixed)
    //! FFTNodes        ( Number of processors for FFT)
    //! Periodic        ( true if periodic BC, false if isolated BC)
    //! Debug           ( If true, lots of output)
    //! GridLength      ( total length of fundamental grid) 
    //!                 ( must be even for periodic simulation)
    //! Padding         (-1) => resolution jumps by factors of 2 (cheap)
    //!                 (1)  => each high point fully padded (expensive)
    //!                 (0)  => no padding (faster, but not recommended)
    //! LevelMax        ( number of levels of refinement (usually 8))
    //! MinimumNumber   ( minimum number of particles belonging to a high point (usually 8))
    //! 
    //! MaxHypreIterations ( max number of iterations in Hypre (20 recommended))
    //! HypreTolerance  ( accuracy of hypre solver (1.0e-7 recommended))
    //! BaseDirectory   ( (string) Base directory for Fractal output, absolute path, ends with /)
    //! RunIdentifier   ( (string) Unique Identifier, no spaces)
    //! TimeTrial       ( If true, lots of MPI_Barrier to check timings)
    //! JMO stuff
    //! G -- the gravitational constant in YOUR units.
    //! xmin -- lower left corner of computational cube in YOUR units.
    //! xmax -- upper right corner of computational cube in YOUR units.
    //! TalkToMe        (Communicator to use for Fractal)
    //
    // Default values. Use YOUR own.
    //
namespace FractalSpace
{
  Fractal_Memory* FractalGravityCreate()
  {
    int NumberParticles=12345;
    int balance=0;
    int FractalNodes0=4;
    int FractalNodes1=4;
    int FractalNodes2=4;
    int FFTNodes=12345;
    bool Periodic=true;
    bool Debug=true;
    int GridLength=256;
    int Padding=-1;
    int LevelMax=8;
    int MinimumNumber=8;
    int MaxHypreIterations=20;
    double HypreTolerance=1.0e-7;
    string BaseDirectory="/p/lscratchc/jensv/";
    string RunIdentifier="SwampCastle";
    bool TimeTrial=true;
    //
    //    double G=6.67428e-11;
    //    vector <double> xmin(3,-2.71828); 
    //    vector <double> xmax(3,1.618034);
    MPI_Comm TalkToMe=MPI_COMM_WORLD;
    //
    Fractal_Memory* PFM=FractalGravityFirstTime(NumberParticles,
						balance,
						FractalNodes0,
						FractalNodes1,
						FractalNodes2,
						FFTNodes,
						Periodic,
						Debug,
						GridLength,
						Padding,
						LevelMax,
						MinimumNumber,
						MaxHypreIterations,
						HypreTolerance,
						BaseDirectory,
						RunIdentifier,
						TimeTrial,
						TalkToMe);
    return PFM;
  }
  void FractalGravity(Fractal_Memory* PFM,int NumberParticles,vector <double>& xmin,vector <double>& xmax, double G)
  {
    PFM->setNumberParticles(NumberParticles);
    fractal_create(PFM);
    //
    int nieach=100; // or whatever.
    int nitotal=PFM->number_particles;
    // "ni" is in honour to the "Knights that say NI".
    vector <double> posx(nieach); // X positions in YOUR units.
    vector <double> posy(nieach); // Y positions in YOUR units.
    vector <double> posz(nieach); // Z positions in YOUR units.
    vector <double> masses(nieach); // masses in YOUR units.
    for(int ni=0;ni<nitotal;ni+=nieach)
      {
	int first=ni;
	int total=min(nieach,nitotal-ni);
	// populate positions and mass vectors with data from YOUR data structure
	// then put the data into *PFM.
	add_particles(PFM,
		      first,
		      total,
		      xmin,
		      xmax,
		      posx,
		      posy,
		      posz,
		      masses);
      }
    do_fractal_force(PFM);
    //
    vector <double>pot; // Potential in YOUR units.
    vector <double>fx;  //  X forces in YOUR units.
    vector <double>fy;  //  Y forces in YOUR units.
    vector <double>fz;  //  Z forces in YOUR units.
    swap(fx,posx);
    swap(fy,posy);
    swap(fz,posz);
    swap(pot,masses);
    //
    for(int ni=0;ni<nitotal;ni+=nieach)
      {
	int first=ni;
	int total=min(nieach,nitotal-ni);
	// Get the potential and forces from *PFM,
	// then put the field into YOUR data structure.
	// G is the gravitational constant in YOUR units.
	// The field is returned in YOUR units.
	// No conversion necessary.
	get_field(PFM,
		  first,
		  total,
		  G,
		  xmin,
		  xmax,
		  pot,
		  fx,
		  fy,
		  fz);
      }
    fractal_delete(PFM);
  }
  void FractalGravityFinal(Fractal_Memory* PFM)
  {
    fractal_memory_content_delete(PFM);
    fractal_memory_delete(PFM);
    PFM=0;
  }
}
