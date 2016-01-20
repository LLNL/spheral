#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
#include "fractal_interface_public.hh"
int main()
{
  using namespace FractalSpace;
  cout << "starting out " << endl;
  Fractal_Memory* PFM=fractal_memory_create();

  int balance=0;
  int NumberParticles=500000;
  int FractalNodes0=1;
  int FractalNodes1=2;
  int FractalNodes2=1;
  int FFTNodes=64;
  bool Periodic=true;
  bool Debug=true;
  int GridLength=256;
  int Padding=-1;
  int LevelMax=8;
  int MinimumNumber=8;
  int MaxHypreIterations=20;
  double HypreTolerance=1.0e-7;
  string BaseDirectory="/p/lscratchd/jensv/";
  string RunIdentifier="TychoBrahe";

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

  fractal_memory_setup(PFM);

  vector <double> posx;
  vector <double> posy;
  vector <double> posz;
  vector <double> velx;
  vector <double> vely;
  vector <double> velz;
  vector <double> masses;

  make_me_a_universe(PFM,posx,posy,posz,velx,vely,velz,masses);

  int FractalNodes=PFM->p_mess->FractalNodes;
  int FractalRank=PFM->p_mess->FractalRank;
  PFM->number_steps_total=50;
  PFM->number_steps_out=100;
  PFM->step_length=1.0e-3;

  for(int step=0;step<PFM->number_steps_total;step++)
    {
      fractal_create(PFM);
      add_particles(PFM,0,NumberParticles,xmin,xmax,posx,posy,posz,masses);
      if(PFM->balance == 1)
	balance_by_particles(PFM);
      do_fractal_force(PFM);
      take_a_leap_cosmo(PFM,masses,G,xmin,xmax,posx,posy,posz,velx,vely,velz);
      am_I_conservative_enough_cosmo(PFM,masses,G,xmin,xmax,-0.5,posx,posy,posz,velx,vely,velz);
      if(step % PFM->number_steps_out == 0)
	start_writing(PFM,NumberParticles,G,xmin,xmax,posx,posy,posz,velx,vely,velz,masses);
      fractal_delete(PFM);
    }
  fractal_memory_content_delete(PFM);
  fractal_memory_delete(PFM);
}
