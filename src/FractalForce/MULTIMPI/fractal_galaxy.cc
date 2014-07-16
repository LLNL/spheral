#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
#include "fractal_interface_public.hh"
int main(int argc, char* argv[])
{
  using namespace FractalSpace;
  MPI_Init(NULL,NULL);
  int FRN;
  MPI_Comm_size(MPI_COMM_WORLD,&FRN);

  bool _inteL_=true;
  if(argc >= 2)
    _inteL_=atoi(argv[1]) != 0;
  int dims[]={0,0,0};
  int GRL=256;
  if(argc >= 3)
    GRL=atoi(argv[2]);
  if(argc >= 4)
    dims[0]=atoi(argv[3]);
  if(argc >= 5)
    dims[1]=atoi(argv[4]);
  if(argc >= 6)
    dims[2]=atoi(argv[5]);
  dims[0]=max(dims[0],0);
  dims[1]=max(dims[1],0);
  dims[2]=max(dims[2],0);
  MPI_Dims_create(FRN,3,dims);
  int FractalNodes0=dims[0];
  int FractalNodes1=dims[1];
  int FractalNodes2=dims[2];
  int NumberParticles=100000;
  if(argc >= 7)
    NumberParticles=atoi(argv[6]);
  cout << "starting out " << argc << " " << FRN << " " << _inteL_ << " " << GRL << " " << FractalNodes0 << " " << FractalNodes1 << " " << FractalNodes2;
  cout << " " << NumberParticles << "\n";

  Fractal_Memory* PFM=fractal_memory_create();
  int balance=1;
  int FFTNodes=9876543;
  bool Periodic=false;
  bool Debug=true;
  int GridLength=GRL;
  if(GRL < 64)
    GridLength=256;
  int Padding=-1;
  int LevelMax=10;
  //  int LevelMax=0;
  int MinimumNumber=8;
  int MaxHypreIterations=20;
  double HypreTolerance=1.0e-7;
  string sa="/p/lscratch";
  string sb="d";
  if(!_inteL_)
    sb="v";
  string sc="/jensv/galaxy/";
  string BaseDirectory=sa+sb+sc;
  //  string BaseDirectory="/p/lscratchc/jensv/galaxy/";
  string RunIdentifier="NerdsRule";
  bool TimeTrial=true;

  FFTNodes=min(FFTNodes,FractalNodes0*FractalNodes1*FractalNodes2);

  PFM->setBalance(balance);
  PFM->setNumberParticles(NumberParticles);
  PFM->setFractalNodes(FractalNodes0,FractalNodes1,FractalNodes2);
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
  fractal_memory_setup(PFM);

  //  int FractalNodes=PFM->p_mess->FractalNodes;
  int FractalRank=PFM->p_mess->FractalRank;
  std::srand(9973+256*FractalRank);
  vector <double> xmin(3,-60.0);
  vector <double> xmax(3,50.0);
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
  PFM->hypre_max_node_load=80000;
  PFM->hypre_max_average_load=20000;
  PFM->number_steps_total=1603;
  //  PFM->number_steps_total=13;
  PFM->number_steps_out=40;
  //  PFM->number_steps_out=200000;
  PFM->step_length=4.0e-5;
  //  PFM->step_length=2.0e-40;
  PFM->time=0.0;
  make_me_a_galaxy(FractalRank,NumberParticles,total_mass,masses,G,posx,posy,posz,velx,vely,velz);
  //  ofstream& FFM=PFM->p_file->FileFractalMemory;
  //  FFM << " info " << NumberParticles << " " << m << " " << total_mass << " " << PFM->time << " " << PFM->step_length << "\n";
  FILE* PFFM=PFM->p_file->PFFractalMemory;
  fprintf(PFFM," info %d %13.4E %13.4E %10.2E %10.2E \n",NumberParticles,m,total_mass,PFM->time,PFM->step_length);
  for(int step=0;step<PFM->number_steps_total;step++)
    {
      fractal_create(PFM);
      add_particles(PFM,0,NumberParticles,xmin,xmax,posx,posy,posz,masses);
      if(PFM->balance > 0)
	balance_by_particles(PFM);
      do_fractal_force(PFM);
      take_a_leap_isol(PFM,masses,G,xmin,xmax,posx,posy,posz,velx,vely,velz);
      am_I_conservative_enough_isol(PFM,masses,G,xmin,xmax,-0.5,posx,posy,posz,velx,vely,velz);
      if(step % PFM->number_steps_out == 0)
	start_writing(PFM,NumberParticles,G,xmin,xmax,posx,posy,posz,velx,vely,velz,masses);
      fractal_delete(PFM);
      //      PFM->p_file->FlushAll();
    }
  fractal_memory_content_delete(PFM);
  fractal_memory_delete(PFM);
}
