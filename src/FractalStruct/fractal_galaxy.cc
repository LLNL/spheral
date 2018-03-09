#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
#include "fractal_interface_public.hh"
int main(int argc, char* argv[])
{
  using namespace FractalSpace;
  int knights;
  MPI_Initialized(&knights);
  if(!knights)
    MPI_Init(NULL,NULL);
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
  double step_length=1.0;
  if(argc >= 11)
    step_length=atof(argv[10]);
  int number_steps_total=1000;
  if(argc >= 12)
    number_steps_total=atoi(argv[11]);
  if(Mess::IAMROOT)
    {
      cerr << "starting out " << argc << " " << FRN << " " << BaseDirectory << " " << RunIdentifier << " " << GRL << " " << FractalNodes0 << " " << FractalNodes1 << " " << FractalNodes2 << "\n";
      cerr << " " << "NumberParticles" << " " << "SHRINK" << " " << "PADDING" << " " << "STEP_LENGTH" << " " << "number_steps_total" << "\n";
      cerr << " " << NumberParticles << " " << SHRINK << " " << PADDING << " " << step_length << " " << number_steps_total << "\n";
      int ar=0;
      while(ar < argc)
	cerr << " " << argv[ar++];
      cerr << "\n";
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
  //  string RunIdentifier="NerdsRule";
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
  PFM->hypre_max_node_load=-1; //////
  PFM->hypre_multiplier=-1;
  fractal_memory_setup(PFM);


  int FractalNodes=PFM->p_mess->FractalNodes;
  FFTNodes=PFM->p_mess->FFTNodes;
  int FractalRank=PFM->p_mess->FractalRank;
  if(FFTNodes < FractalNodes)
    {
      if(PFM->p_mess->IAmAnFFTNode)
	{
	  NumberParticles/=10;
	  PFM->setNumberParticles(NumberParticles);
	}
      PFM->p_mess->calc_total_particles(NumberParticles);
    }
  std::srand(9973+256*FractalRank);
  //  vector <double> xmin(3,-60.0);
  vector <double> xmin(3,-50.0);
  vector <double> xmax(3,50.0);
  vector <double> xmini(3);
  vector <double> xmaxy(3);
  double total_mass=1.0e9;
  double G=3.141592;
  long int TotalNumberParticles=PFM->p_mess->number_particles_total;
  double m=total_mass/static_cast<double>(TotalNumberParticles);
  vector <double> posx(NumberParticles,0.0);
  vector <double> posy(NumberParticles,0.0);
  vector <double> posz(NumberParticles,0.0);
  vector <double> velx(NumberParticles,0.0);
  vector <double> vely(NumberParticles,0.0);
  vector <double> velz(NumberParticles,0.0);
  vector <double> masses(NumberParticles,m);
  PFM->hypre_load_balance=true;
  PFM->number_steps_total=number_steps_total;
  PFM->number_steps_out=20;
  PFM->step_length=step_length;
  PFM->time=0.0;
  make_me_a_galaxy(FractalRank,NumberParticles,total_mass,G,posx,posy,posz,velx,vely,velz);
  FILE* PFFM=PFM->p_file->PFFractalMemory;
  fprintf(PFFM," info %d %d %13.4E %13.4E %10.2E %10.2E \n",TotalNumberParticles,NumberParticles,m,total_mass,PFM->time,PFM->step_length);
  for(int step=0;step<PFM->number_steps_total;step++)
    {
      xmini=xmin;
      xmaxy=xmax;
      shrink_cube(SHRINK,xmin,xmax,PFM,posx,posy,posz,NumberParticles,xmini,xmaxy);
      fractal_create(PFM);
      add_particles(PFM,0,NumberParticles,xmini,xmaxy,posx,posy,posz,masses);
      if(PFM->balance > 0)
	balance_by_particles(PFM,true);
      do_fractal_force(PFM);
      take_a_leap_isol(PFM,masses,G,xmini,xmaxy,posx,posy,posz,velx,vely,velz);
      am_I_conservative_enough_isol(PFM,masses,G,xmini,xmaxy,-0.5,posx,posy,posz,velx,vely,velz);
      if(step % PFM->number_steps_out == 0)
	start_writing(PFM,NumberParticles,G,xmini,xmaxy,posx,posy,posz,velx,vely,velz,masses);
      fractal_delete(PFM);
    }
  fractal_memory_content_delete(PFM);
  fractal_memory_delete(PFM);
  PFM=0;
  return 0;
}
