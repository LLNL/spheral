#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
#include "fractal_interface_public.hh"
namespace FractalSpace
{
  typedef deque<double>::iterator _ITD__;
  Fractal_Memory* FractalGravityFirstTime(
					  bool Periodic,
					  MPI_Comm& TalkToMe,
					  int GridLength,
					  int FractalNodes0,
					  int FractalNodes1,
					  int FractalNodes2,
					  string BaseDirectory,
					  string RunIdentifier
					  )
  {
    Fractal_Memory* PFM=fractal_memory_create();

    PFM->standalone=false;
    PFM->setPeriodic(Periodic);
    PFM->setGridLength(GridLength);
    PFM->setFractalNodes(FractalNodes0,FractalNodes1,FractalNodes2);
    PFM->setBaseDirectory(BaseDirectory);
    PFM->setRunIdentifier(RunIdentifier);
    Fractal_Memory::FRACTAL_UNIVERSE=TalkToMe;
    int Balance=1;
    bool Debug=true;
    int Padding=-1;
    int LevelMax=8;
    int MinimumNumber=8;
    int MaxHypreIterations=20;
    double HypreTolerance=1.0e-7;
    PFM->setBalance(Balance);
    PFM->setDebug(Debug);
    PFM->setPadding(Padding);
    PFM->setLevelMax(LevelMax);
    PFM->setMinimumNumber(MinimumNumber);
    PFM->setHypreIterations(MaxHypreIterations);
    PFM->setHypreTolerance(HypreTolerance);

    return PFM;    
  }
}
namespace FractalSpace
{
  Fractal_Memory* FractalGravityIsolatedFirstTime(
						  MPI_Comm& TalkToMe,
						  int GridLength,
						  int FractalNodes0,
						  int FractalNodes1,
						  int FractalNodes2,
						  string BaseDirectory,
						  string RunIdentifier
						  )
  {
    Fractal_Memory* PFM=fractal_memory_create();

    PFM->standalone=false;
    PFM->setPeriodic(false);
    PFM->setGridLength(GridLength);
    PFM->setFractalNodes(FractalNodes0,FractalNodes1,FractalNodes2);
    PFM->setBaseDirectory(BaseDirectory);
    PFM->setRunIdentifier(RunIdentifier);
    Fractal_Memory::FRACTAL_UNIVERSE=TalkToMe;
    int Balance=1;
    bool Debug=true;
    int Padding=-1;
    int LevelMax=8;
    int MinimumNumber=8;
    int MaxHypreIterations=20;
    double HypreTolerance=1.0e-7;
    PFM->setBalance(Balance);
    PFM->setDebug(Debug);
    PFM->setPadding(Padding);
    PFM->setLevelMax(LevelMax);
    PFM->setMinimumNumber(MinimumNumber);
    PFM->setHypreIterations(MaxHypreIterations);
    PFM->setHypreTolerance(HypreTolerance);

    return PFM;    
  }
  void DoFractalGravity(Fractal_Memory* PFM)
  {
    Fractal* PF=PFM->p_fractal;
    PF->timing(-2,0);
    PF->timing(-1,49);
    fractal_force(*PF,*PFM);
    PFM->steps++;
    PF->timing(1,49);
    PF->timing(0,0);
    PF->timing_lev(0,0);
  }
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
    Fractal_Memory* PFM = new Fractal_Memory;
    return PFM;
  }
  void fractal_memory_delete(Fractal_Memory* PFM)
  {
    if(PFM == 0)
      return;
    delete PFM;
  }
  void fractal_memory_setup(Fractal_Memory* PFM)
  {
    PFM->global_level_max=PFM->level_max;
    /***********************/
    //Leave these parameters alone
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
    PFM->min_hypre_group_size=729;
    //
    PFM->splits=0;
    PFM->masks=0;
    //
    PFM->masks_init=0;
    //
    // Construct a Mess object. 
    // All MPI and FFTW stuff is done in Mess member functions. 
    // This will be used throughout the simulation. 
    Mess* p_mess=new Mess(true,
			  PFM->grid_length,
			  PFM->periodic,
			  PFM->number_particles,
			  PFM->FractalNodes0,
			  PFM->FractalNodes1,
			  PFM->FractalNodes2,
			  PFM->FFTNodes,
			  PFM->FractalWorld);
    PFM->p_mess=p_mess;
    PFM->FFTNodes=PFM->p_mess->FFTNodes;
    p_mess->time_trial=PFM->time_trial;
    p_mess->standalone=PFM->standalone;    
    // Construct a File object. 
    // All output is done in File member functions. 
    // This will be used throughout the simulation.
    File* p_file=0;
    if(PFM->BaseDirectory == "")
      p_file=new File();
    else
      File* p_file=new File(PFM->BaseDirectory,p_mess->FractalNodes,p_mess->FractalRank,PFM->RUN);
    PFM->p_file=p_file;
    PFM->p_mess->p_file=p_file;
    
    // Calculate all simulation information needed. 
    // Includes Boxes, FFTW startup etc.
    // This will be used throughout the simulation.
    PFM->calc_FractalNodes();
    PFM->calc_Buffers_and_more();
    PFM->calc_RealBoxes();
  }
  void fractal_memory_content_delete(Fractal_Memory* PFM)
  {
    delete PFM->p_mess->Parts_in;
    PFM->p_mess->Parts_in=NULL;
    delete PFM->p_fractal;
    PFM->p_fractal=0;
    delete PFM->p_mess;
    PFM->p_mess=0;
    delete PFM->p_file;
    PFM->p_file=0;
  }
  void fractal_create(Fractal_Memory* PFM)
  {
    int NP=PFM->number_particles;
    Fractal* PF=new Fractal(*PFM);
    PFM->p_fractal=PF;
    if(PFM->periodic)
      return;
    Particle* PL;
    try
      {
	PL=new Particle[NP];
      }
    catch(bad_alloc& ba)
      {
	cerr << " bad particles " << NP << " " << ba.what() << endl;
	exit(0);
      }
    PFM->p_mess->Parts_in=PL;
    PF->particle_list.resize(NP);
    for(int ni=0;ni<NP;ni++)
      PF->particle_list[ni]=&PL[ni];
  }
  bool I_am_a_real_particle(Fractal_Memory* PFM,int ni)
  {
    return PFM->p_fractal->particle_list[ni]->get_p_highest_level_group() != 0;
  }
  template <class ForwardIterator>
  void add_particles(Fractal_Memory* PFM,int first,int last,
		     ForwardIterator posxb,ForwardIterator posyb,ForwardIterator poszb,
		     ForwardIterator massesb)
  {
    vector <double> xmin{0.0,0.0,0.0};
    vector <double> xmax{1.0,1.0,1.0};
    add_particles(PFM,first,last,xmin,xmax,posxb,posyb,poszb,massesb);
  }
  template
  void add_particles(Fractal_Memory* PFM,int first,int last,
		     _ITD__ posxb,_ITD__ posyb,_ITD__ poszb,
		     _ITD__ massesb);
  
  template <class ForwardIterator>
  void add_particles(Fractal_Memory* PFM,int first,int last,
		     vector <double> xmin,vector <double> xmax,
		     ForwardIterator posxb,ForwardIterator posyb,ForwardIterator poszb,
		     ForwardIterator massesb)
  {
    vector <double> pos(3);
    double dinv=1.0/(xmax[0]-xmin[0]);
    for(int ni=first;ni<last;ni++)
      {
	pos[0]=(*posxb-xmin[0])*dinv;
	pos[1]=(*posyb-xmin[1])*dinv;
	pos[2]=(*poszb-xmin[2])*dinv;
	PFM->p_fractal->particle_list[ni]->set_posm(pos,*massesb);
	posxb++;
	posyb++;
	poszb++;
	massesb++;
      }
  }
  template
  void add_particles(Fractal_Memory* PFM,int first,int last,
		     vector <double> xmin,vector <double> xmax,
		     _ITD__ posxb,_ITD__ posyb,_ITD__ poszb,
		     _ITD__ massesb);

  void add_particles(Fractal_Memory* PFM,int first,int total,
		     vector <double>& posx,vector <double>& posy,
		     vector <double>& posz,vector <double>& masses)
  {
    vector <double> xmin{0.0,0.0,0.0};
    vector <double> xmax{1.0,1.0,1.0};
    add_particles(PFM,first,total,xmin,xmax,posx,posy,posz,masses);
  }
  void add_particles(Fractal_Memory* PFM,int first,int total,
		     vector <double> xmin,vector <double> xmax,
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
  void return_particles(Fractal_Memory* PFM,int total,
			vector <double> xmin,vector <double> xmax,
			vector <double>& posx,vector <double>& posy,
			vector <double>& posz,vector <double>& masses)
  {
    clean_vector(posx);
    clean_vector(posy);
    clean_vector(posz);
    clean_vector(masses);
    
    double scale=xmax[0]-xmin[0];
    for(int ni=0;ni<total;ni++)
      {
	vector <double> pos(3);
	double mass;
	PFM->p_fractal->particle_list[ni]->get_posm(pos,mass);
	posx.push_back(pos[0]*scale+xmin[0]);
	posy.push_back(pos[1]*scale+xmin[1]);
	posz.push_back(pos[2]*scale+xmin[2]);
	masses.push_back(mass);
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
    unsigned int utotal=total;
    assert(utotal <= pot.size());
    assert(utotal <= fx.size());
    assert(utotal <= fy.size());
    assert(utotal <= fz.size());
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
    delete [] PFM->p_mess->Parts_in;
    PFM->p_mess->Parts_in=0;
    delete PFM->p_fractal;
    PFM->p_fractal=0;
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
  }
  void Fractal_Memory::setFFTNodes(int fmax)
  {
    FFTNodes=fmax;
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
    HTOL=HT;
  }
  void Fractal_Memory::setBaseDirectory(string BD)
  {
    BaseDirectory=BD;
  }
  void Fractal_Memory::setRunIdentifier(string RI)
  {
    RUN=RI;
  }
  void Fractal_Memory::setTimeTrial(bool tt)
  {
    time_trial=tt;
  }
  void Fractal_Memory::setTalkToMe(MPI_Comm& ttm)
  {
    FractalWorld=ttm;
  }
}
