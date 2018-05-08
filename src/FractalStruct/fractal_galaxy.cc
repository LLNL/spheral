#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
int main(int argc, char* argv[])
{
  using namespace FractalSpace;
  int knights;
  MPI_Initialized(&knights);
  if(!knights)
    MPI_Init(NULL,NULL);
  /*
    BaseDirectory
    RunIdentifier
    GridLength
    FractalNodes0
    FractalNodes1
    FractalNodes2
    NumberParticles
    Shrink
    Padding
    NodeLoad
    StepLength
    Number_Steps_Total
  */
  Fractal_Memory::FRACTAL_UNIVERSE=MPI_COMM_WORLD;
  assert(Fractal_Memory::FRACTAL_UNIVERSE != MPI_COMM_NULL);
  int FRN;
  MPI_Comm_size(Fractal_Memory::FRACTAL_UNIVERSE,&FRN);
  int Ranky;
  MPI_Comm_rank(Fractal_Memory::FRACTAL_UNIVERSE,&Ranky);
  Mess::IAMROOT=Ranky == 0;

  string BaseDirectory="/p/lscratchh/jensv/galaxy/";
  if(argc >= 2)
    BaseDirectory=argv[1];

  string RunIdentifier="NerdsRule";
  if(argc >= 3)
    RunIdentifier=argv[2];

  int dims[]={0,0,0};
  int GRL=256;
  if(argc >= 4)
    GRL=atoi(argv[3]);
  if(argc >= 5)
    dims[0]=atoi(argv[4]);
  if(argc >= 6)
    dims[1]=atoi(argv[5]);
  if(argc >= 7)
    dims[2]=atoi(argv[6]);
  dims[0]=max(dims[0],0);
  dims[1]=max(dims[1],0);
  dims[2]=max(dims[2],0);
  MPI_Dims_create(FRN,3,dims);
  int FractalNodes0=dims[0];
  int FractalNodes1=dims[1];
  int FractalNodes2=dims[2];

  int NumberParticles=100000;
  if(argc >= 8)
    NumberParticles=atoi(argv[7]);
  if(NumberParticles < 1000)
    cerr << " TOO FEW PARTICLES " << Ranky << " " << NumberParticles << "\n";

  double SHRINK=0.0;
  if(argc >= 9)
    SHRINK=atof(argv[8]);

  double PADDING=-1;
  if(argc >= 10)
    PADDING=atoi(argv[9]);

  int node_load=20000;
  if(argc >= 11)
    node_load=atoi(argv[10]);

  double step_length=1.0;
  if(argc >= 12)
    step_length=atof(argv[11]);

  int number_steps_total=1000;
  if(argc >= 13)
    number_steps_total=atoi(argv[12]);

  if(Mess::IAMROOT)
    {
      cerr << " Getting Started " << "\n";
      int ar=0;
      while(ar < argc)
	cerr << " " << ar << " " << argv[ar++] << "\n";
    }

  Fractal_Memory* PFM=fractal_memory_create();
  int balance=1;
  int FFTNodes=9876543;
  bool Periodic=false;
  bool Debug=true;
  int GridLength=GRL;
  int LevelMax=10;
  int MinimumNumber=8;
  int MaxHypreIterations=20;
  double HypreTolerance=1.0e-7;
  bool TimeTrial=true;

  FFTNodes=min(FFTNodes,FractalNodes0*FractalNodes1*FractalNodes2);

  PFM->setBalance(balance);
  PFM->setNumberParticles(NumberParticles);
  PFM->setFractalNodes(FractalNodes0,FractalNodes1,FractalNodes2);
  PFM->setFFTNodes(FFTNodes);
  PFM->setPeriodic(Periodic);
  PFM->setDebug(Debug);
  PFM->setGridLength(GridLength);
  PFM->setPadding(PADDING);
  PFM->setLevelMax(LevelMax);
  PFM->setMinimumNumber(MinimumNumber);
  PFM->setHypreIterations(MaxHypreIterations);
  PFM->setHypreTolerance(HypreTolerance);
  PFM->setBaseDirectory(BaseDirectory);
  PFM->setRunIdentifier(RunIdentifier);
  PFM->setTimeTrial(TimeTrial);
  PFM->hypre_max_node_load=node_load;
  PFM->hypre_multiplier=-1;
  fractal_memory_setup(PFM);


  int FractalNodes=PFM->p_mess->FractalNodes;
  FFTNodes=PFM->p_mess->FFTNodes;
  int FractalRank=PFM->p_mess->FractalRank;
  if(FFTNodes < FractalNodes)
    {
      if(PFM->p_mess->IAmAnFFTNode)
	{
	  // NumberParticles/=10;
	  PFM->setNumberParticles(NumberParticles);
	}
      PFM->p_mess->calc_total_particles(NumberParticles);
    }
  std::srand(9973+256*FractalRank);
  vector <double> xmin(3,-50.0);
  vector <double> xmax(3,50.0);
  vector <double> xmini(3);
  vector <double> xmaxy(3);
  double total_mass=1.0e9;
  double G=3.141592;
  long int TotalNumberParticles=PFM->p_mess->number_particles_total;
  double m=total_mass/static_cast<double>(TotalNumberParticles);
  deque <double> posx(NumberParticles,0.0);
  deque <double> posy(NumberParticles,0.0);
  deque <double> posz(NumberParticles,0.0);
  deque <double> velx(NumberParticles,0.0);
  deque <double> vely(NumberParticles,0.0);
  deque <double> velz(NumberParticles,0.0);
  deque <double> masses(NumberParticles,m);
  PFM->hypre_load_balance=true;
  PFM->number_steps_total=number_steps_total;
  PFM->number_steps_out=20;
  PFM->step_length=step_length;
  PFM->time=0.0;
  
  make_me_a_galaxy(PFM,
		   NumberParticles,
		   total_mass,
		   G,
		   posx.begin(),
		   posy.begin(),
		   posz.begin(),
		   velx.begin(),
		   vely.begin(),
		   velz.begin());
  FILE* PFFM=PFM->p_file->PFFractalMemory;
  fprintf(PFFM," info %d %d %13.4E %13.4E %10.2E %10.2E \n",TotalNumberParticles,NumberParticles,m,total_mass,PFM->time,PFM->step_length);
  for(int step=0;step<PFM->number_steps_total;step++)
    {
      xmini=xmin;
      xmaxy=xmax;

      shrink_cube(SHRINK,
		  xmin,
		  xmax,
		  PFM,
		  posx.begin(),
		  posy.begin(),
		  posz.begin(),
		  NumberParticles,
		  xmini,
		  xmaxy);

      fractal_create(PFM);

      add_particles(PFM,
		    0,NumberParticles,
		    xmini,
		    xmaxy,
		    posx.begin(),
		    posy.begin(),
		    posz.begin(),
		    masses.begin());
      
      if(PFM->balance > 0)
	balance_by_particles(PFM,true);

      do_fractal_force(PFM);

      take_a_leap_isol(PFM,
		       masses.begin(),
		       G,
		       xmini,
		       xmaxy,
		       posx.begin(),
		       posy.begin(),
		       posz.begin(),
		       velx.begin(),
		       vely.begin(),
		       velz.begin());
      
      am_I_conservative_enough_isol(PFM,
				    masses.begin(),
				    G,
				    xmini,
				    xmaxy,
				    -0.5,
				    posx.begin(),
				    posy.begin(),
				    posz.begin(),
				    velx.begin(),
				    vely.begin(),
				    velz.begin());
      
      if(step % PFM->number_steps_out == 0)
	start_writing(PFM,
		      NumberParticles,
		      G,
		      xmini,
		      xmaxy,
		      posx.begin(),
		      posy.begin(),
		      posz.begin(),
		      velx.begin(),
		      vely.begin(),
		      velz.begin(),
		      masses.begin());
      
      fractal_delete(PFM);
    }
  fractal_memory_content_delete(PFM);
  
  fractal_memory_delete(PFM);
  
  PFM=0;
  return 0;
}
