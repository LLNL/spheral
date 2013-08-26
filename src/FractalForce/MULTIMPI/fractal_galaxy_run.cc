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
  int NumberParticles=500000;
  int FractalNodes0=4;
  int FractalNodes1=4;
  int FractalNodes2=4;
  int FFTNodes=9876543;
  bool Periodic=false;
  bool Debug=true;
  int GridLength=256;
  int Padding=-1;
  int LevelMax=8;
  int MinimumNumber=8;
  int MaxHypreIterations=20;
  double HypreTolerance=1.0e-7;
  string BaseDirectory="/p/lscratchc/jensv/";
  string RunIdentifier="KongenErEnFinke";
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

  int FractalNodes=PFM->p_mess->FractalNodes;
  int FractalRank=PFM->p_mess->FractalRank;
  std::srand(9973+256*FractalRank);
  vector <double> xmin(3,-60.0);
  vector <double> xmax(3,50.0);
  double total_mass=1.0e9;
  double G=3.141592;
  int TotalNumberParticles=PFM->p_mess->number_particles_total;
  double m=total_mass/static_cast<double>(TotalNumberParticles);
  vector <double> posx(NumberParticles,0.0);
  vector <double> posy(NumberParticles,0.0);
  vector <double> posz(NumberParticles,0.0);
  vector <double> velx(NumberParticles,0.0);
  vector <double> vely(NumberParticles,0.0);
  vector <double> velz(NumberParticles,0.0);
  vector <double> masses(NumberParticles,m);
  PFM->number_steps_total=500;
  PFM->number_steps_out=20;
  PFM->step_length=1.0e-5;
  PFM->time=0.0;
  make_me_a_galaxy(NumberParticles,total_mass,masses,G,posx,posy,posz,velx,vely,velz);
  ofstream& FFM=PFM->p_file->FileFractalMemory;
  FFM << " info " << NumberParticles << " " << m << " " << total_mass << " " << PFM->time << " " << PFM->step_length << endl;
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
    }
  fractal_memory_content_delete(PFM);
  fractal_memory_delete(PFM);
}
