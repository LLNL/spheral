#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
#include "fractal_interface_public.hh"
namespace FractalSpace
{
  //! sequence of calls

  //! Setup
  //! Fractal_Memory* PFM=fractal_memory_create();
  //!
  //! PFM->setBalance;
  //! PFM->setNumberParticles;
  //! PFM->setFractalNodes;
  //! PFM->setPeriodic;
  //! PFM->setDebug;
  //! PFM->setGridLength;
  //! PFM->setPadding;
  //! PFM->setLevelMax;
  //! PFM->setMinimumNumber;
  //! PFM->setHypreIterations;
  //! PFM->setHypreTolerance;
  //! PFM->setBaseDirectory;
  //! PFM->setRunIdentifier;
  //!
  //! PFM->fractal_memory_setup();
  //!
  //! for(int ni=0;ni<steps;ni++)
  //! {
  //! Do Your own stuff;
  //! fractal_create;
  //! addParticles;
  //! doFractalForce;
  //! getField;
  //! fractal_delete;
  //! Do Your own stuff;
  //! }
  //! Finish up.
  //! fractal_memory_content_delete;
  //! fractal_memory_delete;



  void do_fractal_force(Fractal_Memory* PFM)
  {
    Fractal* PF=PFM->p_fractal;
    PF->set_steps(PFM->steps);
    PF->timing(-2,0);
    PF->timing(-1,49);
    fractal_force(*PF,*PFM);
    PFM->steps++;
    PF->timing(1,49);
    PF->timing(0,0);
    PF->timing_lev(0,0);
  }
  Fractal_Memory* fractal_memory_create()
  {
    Fractal_Memory* PFM= new Fractal_Memory;
    return PFM;
  }
  void fractal_memory_delete(Fractal_Memory* PFM)
  {
    delete PFM;
  }
  void fractal_memory_setup(Fractal_Memory* PFM)
  {
    PFM->MPIrun=PFM->FractalNodes > 1;
    PFM->global_level_max=PFM->level_max;
    /***********************/
    //Leave these parameters alone
    PFM->min_hypre_group_size=1;
    PFM->new_points_gen=9;
    PFM->steps=0;
    PFM->momentum_conserve=false;
    PFM->amnesia=true; // (true) forget everything after you are done. (false) remember everything.
    PFM->mind_wipe=false; // (true) delete everything and then come back without calculating anything.
    PFM->fixed_potential=false; // (true) use the fixed potential.
    PFM->calc_shear=false;// (true) if we calculate shear of force field
    PFM->calc_density_particle=false;
    PFM->do_vel=false;
    PFM->start_up=false;
    PFM->halo_fixed=false;
    //
    PFM->splits=0;
    PFM->masks=0;
    //
    PFM->masks_init=0;
    //
    // Construct a Mess object. 
    // All MPI and FFTW stuff is done in Mess member functions. 
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
    PFM->p_file->note(true," a fractal_memory ");
    PFM->calc_FractalNodes();
    PFM->p_file->note(true," b fractal_memory ");
    PFM->calc_Buffers_and_more();
    PFM->p_file->note(true," c fractal_memory ");
    PFM->calc_RealBoxes();
    PFM->p_file->note(true," d fractal_memory ");
  }
  void fractal_memory_content_delete(Fractal_Memory* PFM)
  {
    Fractal* PF=PFM->p_fractal;
    Particle* P=PFM->p_mess->Parts_in;
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
  void fractal_create(Fractal_Memory* PFM)
  {
    int NP=PFM->number_particles;
    Fractal* PF=new Fractal(*PFM);
    PFM->p_fractal=PF;
    Particle* PL=new Particle[NP];
    PFM->p_mess->Parts_in=PL;
    PF->particle_list.resize(NP);
    for(int ni=0;ni<NP;ni++)
      PF->particle_list[ni]=&PL[ni];
  }
  bool I_am_a_real_particle(Fractal_Memory* PFM,int ni)
  {
    return PFM->p_fractal->particle_list[ni]->get_p_highest_level_group() != 0;
  }
  void add_particles(Fractal_Memory* PFM,int first,int total,
		     vector <double>& xmin,vector <double>& xmax,
		     vector <double>& posx,vector <double>& posy,
		     vector <double>& posz,vector <double>& masses)
  {
    total=min(first+total,PFM->number_particles)-first;
    vector <double> pos(3);
    double dinv=1.0/(xmax[0]-xmin[0]);
    for(int ni=0;ni<total;ni++)
      {
	pos[0]=(posx[ni]-xmin[0])*dinv;
	pos[1]=(posy[ni]-xmin[1])*dinv;
	pos[2]=(posz[ni]-xmin[2])*dinv;
	PFM->p_fractal->particle_list[ni+first]->set_posm(pos,masses[ni]);
      }
  }
  void get_potential(Fractal_Memory* PFM,int first,int total,double G,
		     vector <double>& xmin,vector <double>& xmax,vector <double>& pot)
  {
    total=min(first+total,PFM->number_particles)-first;
    double convpot=G/(xmax[0]-xmin[0]);
    for(int ni=0;ni<total;ni++)
      pot[ni]=PFM->p_fractal->particle_list[ni+first]->get_potential()*convpot;
  }
  void get_field(Fractal_Memory* PFM,int first,int total,double G,
		vector <double>& xmin,vector <double>& xmax,
		vector <double>& pot,vector <double>& fx,
		vector <double>& fy,vector <double>& fz)
  {
    total=min(first+total,PFM->number_particles)-first;
    assert(total <= pot.size());
    assert(total <= fx.size());
    assert(total <= fy.size());
    assert(total <= fz.size());
    vector <double> potforce(4);
    double dinv=1.0/(xmax[0]-xmin[0]);
    double convpot=G*dinv;
    double convforce=G*dinv*dinv;
    for(int ni=0;ni<total;ni++)
      {
	PFM->p_fractal->particle_list[ni+first]->get_field_pf(potforce);
	pot[ni]=potforce[0]*convpot;
	fx[ni]=potforce[1]*convforce;
	fy[ni]=potforce[2]*convforce;
	fz[ni]=potforce[3]*convforce;
      }
  }
  void fractal_delete(Fractal_Memory* PFM)
  {
    Fractal* PF=PFM->p_fractal;
    Particle* P=PFM->p_mess->Parts_in;
    delete [] P;
    P=0;
    delete PF;
    PF=0;
  }
  void Fractal_Memory::setBalance(int B)
  {
    balance=B;
  }
  void Fractal_Memory::setNumberParticles(int NP)
  {
    number_particles=NP;
  }
  void Fractal_Memory::setFractalNodes(int FR0,int FR1,int FR2)
  {
    FractalNodes0=FR0;
    FractalNodes1=FR1;
    FractalNodes2=FR2;
    FractalNodes=FR0*FR1*FR2;
    MPIrun=FractalNodes > 1;
  }
  void Fractal_Memory::setPeriodic(bool per)
  {
    periodic=per;
  }
  void Fractal_Memory::setDebug(bool Db)
  {
    debug=Db;
  }
  void Fractal_Memory::setGridLength(int GL)
  {
    grid_length=GL;
  }
  void Fractal_Memory::setPadding(int PA)
  {
    PA=min(PA,1);
    PA=max(PA,-1);
    padding=PA;
  }
  void Fractal_Memory::setLevelMax(int LM)
  {
    level_max=LM;
  }
  void Fractal_Memory::setMinimumNumber(int MN)
  {
    minimum_number=MN;
  }
  void Fractal_Memory::setHypreIterations(int MHI)
  {
    maxits=MHI;
  }
  void Fractal_Memory::setHypreTolerance(double HT)
  {
    epsilon_sor=HT;
  }
  void Fractal_Memory::setBaseDirectory(string BD)
  {
    BaseDirectory=BD;
  }
  void Fractal_Memory::setRunIdentifier(string RI)
  {
    RUN=RI;
  }
}
