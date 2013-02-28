#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
#include "fractal_interface_public.hh"
int main()
{
  using namespace FractalSpace;
  cout << "starting out " << endl;
  Fractal_Memory* PFM=fractal_memory_create();

  int balance=1;
  int NumberParticles=200000;
  int FractalNodes0=2;
  int FractalNodes1=3;
  int FractalNodes2=4;
  bool Periodic=false;
  bool Debug=true;
  int GridLength=128;
  int Padding=-1;
  int LevelMax=8;
  int MinimumNumber=8;
  int MaxHypreIterations=20;
  double HypreTolerance=1.0e-7;
  string BaseDirectory="/p/lscratchc/jensv/";
  string RunIdentifier="KongenErEnFinke";

  PFM->setBalance(balance);
  PFM->setNumberParticles(NumberParticles);
  PFM->setFractalNodes(FractalNodes0,FractalNodes1,FractalNodes2);
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

  fractal_memory_setup(PFM);

  int FractalNodes=PFM->p_mess->FractalNodes;
  int FractalRank=PFM->p_mess->FractalRank;
  std::srand(9973+256*FractalRank);
  vector <double> xmin(3,-50.0);
  vector <double> xmax(3,50.0);
  double total_mass=1.0e7;
  double G=2.718281828;
  double m=total_mass/static_cast<double>(NumberParticles*FractalNodes);
  vector <double> posx(NumberParticles,0.0);
  vector <double> posy(NumberParticles,0.0);
  vector <double> posz(NumberParticles,0.0);
  vector <double> velx(NumberParticles,0.0);
  vector <double> vely(NumberParticles,0.0);
  vector <double> velz(NumberParticles,0.0);
  vector <double> masses(NumberParticles,m);
  PFM->number_steps_total=2003;
  PFM->number_steps_out=100;
  PFM->step_length=1.0e-4;
  PFM->time=0.0;
  make_me_a_galaxy(NumberParticles,total_mass,masses,G,posx,posy,posz,velx,vely,velz);

  for(int step=0;step<PFM->number_steps_total;step++)
    {
      fractal_create(PFM);
      add_particles(PFM,0,NumberParticles,xmin,xmax,posx,posy,posz,masses);
      if(PFM->balance == 1)
	balance_by_particles(PFM);
      do_fractal_force(PFM);
      take_a_leap_isol(PFM,masses,G,xmin,xmax,posx,posy,posz,velx,vely,velz);
      am_I_conservative_enough_isol(PFM,masses,G,xmin,xmax,-0.5,posx,posy,posz,velx,vely,velz);
      if(step % PFM->number_steps_out == 0)
	start_writing(PFM,NumberParticles,G,xmin,xmax,posx,posy,posz,velx,vely,velz,masses);
      fractal_delete(PFM);
    }
  fractal_memory_content_delete(PFM);
  fractal_memory_delete(PFM);
}
