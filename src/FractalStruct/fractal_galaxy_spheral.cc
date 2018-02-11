#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
#include "fractal_interface_public.hh"
int main(int argc, char* argv[])
{
  int knights;
  MPI_Initialized(&knights);
  if(!knights)
    MPI_Init(NULL,NULL);
  MPI_Comm& TalkToMe=MPI_COMM_WORLD;
  int GridLength=256;
  if(argc >= 2)
    GridLength=atoi(argv[1]);
  int FractalNodes0=8;
  if(argc >= 3)
    FractalNodes0=atoi(argv[2]);
  int FractalNodes1=8;
  if(argc >= 4)
    FractalNodes0=atoi(argv[3]);
  int FractalNodes1=8;
  if(argc >= 5)
    FractalNodes2=atoi(argv[4]);
  int FRN=-1;
  MPI_Comm_size(MPI_COMM_WORLD,&FRN);
  int RANK=-1;
  MPI_Comm_rank(MPI_COMM_WORLD,&RANK);
  assert(FractalNodes0*FractalNodes1*FractalNodes2 == FRN);
  string Basedirectory="/p/lscratchh/jensv/galaxy/";
  if(argc >= 6)
    Basedirectory=argv[5];
  string RunIdentifier="NerdsRule";
  if(argc >= 7)
    RunIdentifier=argv[6];
  int NumberParticles=100000;
  if(argc >= 8)
    NumberParticles=atoi(argv[7]);
  double G=123.4;
  if(argc >= 9)
    G=atof(argv[8]);
  double slope=-1.5;
  if(argc >= 10)
    slope=atof(argv[8]);
  double RMAX=76.5;
  if(argc >= 10)
    RMAX=atof(argv[9]);
  double x0=16.5;
  if(argc >= 11)
    x0=atof(argv[10]);
  double y0=17.5;
  if(argc >= 12)
    y0=atof(argv[11]);
  double z0=15.5;
  if(argc >= 13)
    z0=atof(argv[12]);
  double mm=15.5;
  if(argc >= 14)
    mm=atof(argv[13]);
  double BOXLENGTH=2.2*RMAX;
  if(argc >= 15)
    BOXLENGTH=atof(argv[14]);
  double SHRINK=1.0;
  if(argc >= 16)
    SHRINK=atof(argv[15]);
  int RANDOMSEED=8765;
  if(argc >=17)
    RANDOMSEED=atoi(argv[16]);
  int withparts=1;
  if(argc >=18)
    withparts=atoi(argv[17]);
  
  srand(withparts+RANK);
  double RMAX2=RMAX/2;
  vector<double>xmin{x0-RMAX2,y0-RMAX2,z0-RMAX2);
  vector<double>xmax{x0+RMAX2,y0+RMAX2,z0+RMAX2);
  vector<double>xmini=xmin;
  vector<double>xmaxy=xmax;
  
  FractalSpace::Fractal_Memory* PFM=
    FractalSpace::FractalGravityIsolatedFirstTime(
						  MPI_Comm& TalkToMe,
						  int GridLength,
						  int FractalNodes0,
						  int FractalNodes1,
						  int FractalNodes2,
						  string BaseDirectory,
						  string RunIdentifier);

  FractalSpace::fractal_memory_setup(FractalSpace::PFM);
  FractalSpace::PFM->setNumberParticles(NumberParticles);
  FractalSpace::fractal_create(FractalSpace::PFM);

  double slope3=slope+3.0;
  double expo=1.0/slope3;
  double rand_max=(double)RAND_MAX;
  double twopi=8.0*atan(1.0);
  vector<double>posx;
  vector<double>posy;
  vector<double>posz;
  vector<double>accx0;
  vector<double>accy0;
  vector<double>accz0;
  double m0=mm/(double)(FRN*NumberParticles);
  for(int ni=0;ni<NumberParticles;ni++)
    {
      double r0=FractalSpace::Fractal::my_rand_not_zero(rand_max);
      double r1=RMAX*pow(r0,expo);
      double zr=2.0*Fractal::my_rand(rand_max)-1.0;
      double Rr=sqrt(max(1.0-zr*zr,0.0));
      double phi=twopi*FractalSpace::Fractal::my_rand(rand_max);
      double xr=Rr*cos(phi);
      double yr=Rr*sin(phi);
      posx.push_back(r1*xr+x0);
      posy.push_back(r1*yr+y0);
      posz.push_back(r1*zr+z0);
      masses.push_back(m0);
      double massr=mm*pow(r0,slope3);
      double Gr2minv=G*massr/(r1*r1);
      accx0.push_back(-Gr2minv*xr);
      accy0.push_back(-Gr2minv*yr);
      accz0.push_back(-Gr2minv*zr);
    }

  FractalSpace::add_particles(FractalSpace::PFM,0,NumberParticles,posx,posy,posz,masses);
  
  FractalSpace::FractalCube(FractalSpace::PFM,SHRINK,xmin,xmax,xmini,xmaxy);

  FractalSpace::balance_by_particles(FractalSpace::PFM,withparts);

  FractalSpace::DoFractalGravity(FractalSpace::PFM);

  vector<double>pot;
  vector<double>accx;
  vector<double>accy;
  vector<double>accz;
  FractalSpace::get_field(FractalSpace::PFM,0,NumberParticles,G,
			  vector <double>& xmini,vector <double>& xmaxy,
			  vector <double>& pot,vector <double>& accx,
			  vector <double>& accy,vector <double>& accz);

  FractalSpace::FILE* PFPos=FractalSpace::PFM->p_file->PFPos;
  for(int ni=0;ni<NumberParticles;ni++)
    {
      FilePos << scientific << ni << "\t" << accx[ni] << "\t" << accx0[ni] << "\t";
      FilePos << scientific << accy[ni] << "\t" << accy0[ni];
      FilePos << scientific << accz[ni] << "\t" << accz0[ni] << "\n";
    }
  
  FractalSpace::fractal_delete(FractalSpace::PFM);  
  
  FractalSpace::fractal_memory_content_delete(FractalSpace::PFM);

  FractalSpace::fractal_memory_delete(FractalSpace::PFM);

  FractalSpace::PFM=0;
  
  return 0;
}
