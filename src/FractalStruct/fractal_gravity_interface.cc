
// Include Spheral include files
#include "whatever_spheral.hh"

// Include Fractal include files
#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
#include "fractal_interface_public.hh"

namespace Spheral 
{
  namespace GravitySpace 
  {
    using namespace FractalSpace;
    using YOUR_SPACES::WhatEver;
    template <T> void fractal_gravity_interface(T& JMO)
    {
      static bool done_fractal_memory=false;
      Fractal_Memory* PFM=0;
      if(JMO.it_is_all_over)
	{
	  FractalSpace::fractal_memory_content_delete(PFM);
	  FractalSpace::fractal_memory_delete(PFM);
	  done_fractal_memory=false;
	  PFM=0;
	  return;
	}
      if(!done_fractal_memory)
	{
	  JMO.p_fractal_memory=FractalSpace::fractal_memory_create();
	  PFM=JMO.p_fractal_memory;
	  //
	  // You have to set these parameters
	  // The Gravitational constant
	  double G=1.0;

	  // Lower left corner of the computational CUBE
	  vector <double>xmin(3,0.0);

	  // Upper right corner of the computational CUBE
	  // Must be a cube
	  vector <double>xmax(3,1.0);

	  // Nodes in each dimension. You can set this or just leave this alone and let
	  // the code figure it out.
	  vector <int>NODEdims(3,0);

	  // Loadbalance nodes.
	  // 0 means equal volume boxes.
	  // 1 means equal number of particles in each box.
	  int balance=1;

	  // Number of particles stored on this particular node.
	  int NumberParticles=100000;

	  // Length of Fourier Grid. Must be even, should have only primefactors <= 5.
	  int GridLength=512;

	  // What is your communicator?
	  MPI_Comm SpheralWorld=MPI_COMM_WORLD;

	  //
	  // You need to specify where to put output files.
	  // In this example Task "N" will put all files in /p/lscratchc/jensv/galaxy/NerdsRuleA_N/;
	  // If "A" is already used, program will use "B" and so on.
	  string BaseDirectory="/p/lscratchc/jensv/galaxy/";
	  string RunIdentifier="NerdsRule";
	  //
	  // Leave these parameters alone for now
	  int FFTNodes=9876543;
	  bool Periodic=false;
	  bool Debug=true;
	  int Padding=-1;
	  int LevelMax=10;
	  int MinimumNumber=8;
	  int MaxHypreIterations=20;
	  double HypreTolerance=1.0e-7;
	  bool TimeTrial=true;
	  //
	  PFM->set_G(G);
	  PFM->set_xmin(xmin);
	  PFM->set_xmax(xmax);
	  PFM->setBalance(balance);
	  PFM->setNumberParticles(NumberParticles);
	  PFM->setFFTNodes(FFTNodes);
	  PFM->setPeriodic(Periodic);
	  PFM->setDebug(Debug);
	  PFM->setGridLength(GridLength);
	  PFM->setPadding(Padding);
	  PFM->setLevelMax(LevelMax);
	  PFM->setMinimumNumber(MinimumNumber);
	  PFM->setHypreIterations(MaxHypreIterations);
	  PFM->setHypreTolerance(HypreTolerance);
	  PFM->setBaseDirectory(BaseDirectory);
	  PFM->setRunIdentifier(RunIdentifier);
	  PFM->setTimeTrial(TimeTrial);
	  FractalSpace::fractal_memory_setup(PFM);
	  PFM->p_mess->createFractalWorld(SpheralWorld,NODEdims);
	  PFM->setFractalNodes(NODEdims[0],NODEdims[1],NODEdims[2]);
	  FFTNodes=min(FFTNodes,NODEdims[0]*NODEdims[1]*NODEdims[2]);
	  done_fractal_memory=true;
	}
      PFM=JMO.p_fractal_memory;
      fractal_create(PFM);
      PFM->get_G(G);
      PFM->get_xmin(xmin);
      PFM->get_xmax(xmax);
      int nump=1000;
      vector <double> posx(nump);
      vector <double> posy(nump);
      vector <double> posz(nump);
      vector <double> masses(nump);
      int number0=0;
      int number1=number0+nump-1;
      while(number0 < NumberParticles)
	{
	  int number1a=min(number1,NumberParticles);

	  // You need to write some code to extract particle positions and masses from your data structures.
	  // You do not need to do any data conversion.
	  get_particle_position_masses_from_spheral(JMO,number0,number1a,posx,posy,posz,masses);
	  // code will stick this data into my data structures

	  FractalSpace::add_particles(PFM,number0,number1a,xmin,xmax,posx,posy,posz,masses);
	  number0+=nump;
	  number1+=nump;
	}
      if(PFM->balance > 0)
	FractalSpace::balance_by_particles(PFM);
      FractalSpace::do_fractal_force(PFM);
      vector <double> pot(nump);
      vector <double> fx(nump);
      vector <double> fy(nump);
      vector <double> fy(nump);
      number0=0;
      number1=number0+nump-1;
      while(number0 < NumberParticles)
	{
	  int number1a=min(number1,NumberParticles);
	  FractalSpace::get_field(PFM,number0,number1a,NumberParticles,
				  G,xmin,xmax,pot,fx,fy,fz);

	  // You need to write some code to stick the potential and accelation back into your data structures.
	  // Data is returned to you in YOUR units.

	  put_particle_field_to_spheral(JMO,number0,number1a,xmin,xmax,pot,fx,fy,fz);
	  number0+=nump;
	  number1+=nump;
	}
      FractalSpace::fractal_delete(PFM);
    }
  }
}
