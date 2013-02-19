#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
#include "fractal_interface_public.hh"
int main()
{
  using namespace FractalSpace;
  cout << "starting out " << endl;
  Fractal_Memory* PFM=fractal_memory_create();

  int NumberParticles=1000000;
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
  double HypreIterations=1.0e-7;
  string BaseDirectory="/p/lscratchc/jensv/";
  string RunIdentifier="KongenErEnFinke";

  PFM->setNumberParticles(NumberParticles);
  PFM->setFractalNodes(FractalNodes0,FractalNodes1,FractalNodes2);
  PFM->setPeriodic(Periodic);
  PFM->setDebug(Debug);
  PFM->setGridLength(GridLength);
  PFM->setPadding(Padding);
  PFM->setLevelMax(LevelMax);
  PFM->setMinimumNumber(MinimumNumber);
  PFM->setHypreIterations(MaxHypreIterations);
  PFM->setHypreTolerance(HypreIterations);
  PFM->setBaseDirectory(BaseDirectory);
  PFM->setRunIdentifier(RunIdentifier);

  fractal_memory_setup(PFM);

  int FractalNodes=PFM->p_mess->FractalNodes;
  int FractalRank=PFM->p_mess->FractalRank;
  int RandomNumber=9973+256*FractalRank;
  std::srand(RandomNumber);
  vector <double> xmin(3,-50.0);
  vector <double> xmax(3,50.0);
  double total_mass=1.0e7;

  double m=total_mass/static_cast<double>(NumberParticles*FractalNodes);
  vector <double> posx(NumberParticles,0.0);
  vector <double> posy(NumberParticles,0.0);
  vector <double> posz(NumberParticles,0.0);
  vector <double> velx(NumberParticles,0.0);
  vector <double> vely(NumberParticles,0.0);
  vector <double> velz(NumberParticles,0.0);

  MakeMeaGalaxy(NumberParticles,posx,posy,posy,velx,vely,velz);

  vector <double>masses(NumberParticles,m);
  for(int step=0;step<PFM->number_steps_total;step++)
    {
      fractal_create(PFM);
      addParticles(PFM,0,NumberParticles,xmin,xmax,posx,posy,posy,masses);
      doFractalForce(PFM);
      takeALeapIsol(PFM,masses,posx,posy,posz,velx,vely,velz);
      fractal_delete(PFM);
    }
  fractal_memory_content_delete(PFM);
  fractal_memory_delete(PFM);
}
#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void takeALeapIsol(Fractal_Memory* PFM,vector <double>& masses,
		     vector <double>& xmin,vector <double>& xmax,
		     vector <double>& posx,vector <double>& posy,vector <double>& posz,
		     vector <double>& velx,vector <double>& vely,vector <double>& velz)
  {
    int stride=100;
    int NP=PFM->number_particles;
    int dt=PFM->step_length;
    vector <double>pot(stride);
    vector <double>fx(stride);
    vector <double>fy(stride);
    vector <double>fz(stride);
    double pe=0.0;
    double ke=0.0;
    double px=0.0;
    double py=0.0;
    double pz=0.0;
    for(int ni=0;ni<NP;ni+=stride)
      {
	int many=min(ni+stride,NP)-ni;
	getField(PFM,ni,ni+many,1.0,xmin,xmax,pot,fx,fy,fz);
	for(int p=0;p<many;p++)
	  {
	    int nip=ni+p;
	    velx[nip]+=fx[p]*dt;
	    vely[nip]+=fy[p]*dt;
	    velz[nip]+=fz[p]*dt;
	    posx[nip]+=velx[nip]*dt;
	    posy[nip]+=vely[nip]*dt;
	    posz[nip]+=velz[nip]*dt;
	    pe+=masses[nip]*pot[p];
	    ke+=masses[nip]*(velx[nip]*velx[nip]+vely[nip]*vely[nip]+velz[nip]*velz[nip]);
	    px+=masses[nip]*velx[nip];
	    py+=masses[nip]*vely[nip];
	    pz+=masses[nip]*velz[nip];
	  }
      }
  }
}
#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void MakeMeaGalaxy(int numbers,
		     vector <double>& posx,vector <double>& posy,vector <double>& posz,
		     vector <double>& velx,vector <double>& vely,vector <double>& velz)
  {
    double rand_max=(double)RAND_MAX;
    double rmax=30.0;
    double x_off=10.0;
    double y_off=-13.0;
    double z_off=11.0;
    double slope=-0.5;
    double total_mass=1.0e7;
    double velratio=0.5;
    double slope3=slope+3.0;
    double expo=1.0/(3.0+slope);
    double twopi=8.0*atan(1.0);
    for(int ni=0;ni<numbers;ni++)
      {
	double r1=Fractal::my_rand(rand_max);
	r1=pow(r1,expo);
	double ctheta=2.0*Fractal::my_rand(rand_max)-1.0;
	double phi=twopi*Fractal::my_rand(rand_max);
	double r=rmax*r1;
	posz[ni]=r*ctheta+z_off;
	double stheta=sqrt(1.0-ctheta*ctheta);
	double cphi=cos(phi);
	double sphi=sin(phi);
	posx[ni]=r*stheta*cphi+x_off;
	posy[ni]=r*stheta*sphi+y_off;
	double massr=total_mass*pow((r/rmax),slope3);
	double vt=velratio*sqrt(massr/r);
	velx[ni]=-sphi*vt;
	vely[ni]=cphi*vt;
	velz[ni]=0.0;
      }
  }
}
