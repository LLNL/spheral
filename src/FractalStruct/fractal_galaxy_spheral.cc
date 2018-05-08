#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
#include "fractal_interface_public.hh"
int main(int argc, char* argv[])
{
  /*
    int GridLength
    int FractaNodes0
    int FractaNodes1
    int FractaNodes2
    string BaseDirectory
    string RunIdentifier
    int NumberParticles
    double G
    double slope
    double RMAX
    double x0
    double y0
    double z0
    double mm
    double BOXLENGTH
    double SHRINK
    int RANDOMSEED
    int withparts
  */
  int knights;
  MPI_Initialized(&knights);
  if(!knights)
    MPI_Init(NULL,NULL);
  MPI_Comm TalkToMe=MPI_COMM_WORLD;
  int GridLength=256;
  if(argc >= 2)
    GridLength=atoi(argv[1]);
  int FractalNodes0=8;
  if(argc >= 3)
    FractalNodes0=atoi(argv[2]);
  int FractalNodes1=8;
  if(argc >= 4)
    FractalNodes1=atoi(argv[3]);
  int FractalNodes2=8;
  if(argc >= 5)
    FractalNodes2=atoi(argv[4]);
  int FRN=-1;
  MPI_Comm_size(MPI_COMM_WORLD,&FRN);
  int RANK=-1;
  MPI_Comm_rank(MPI_COMM_WORLD,&RANK);
  assert(FractalNodes0*FractalNodes1*FractalNodes2 == FRN);
  string BaseDirectory="/p/lscratchh/jensv/galaxy/";
  if(argc >= 6)
    BaseDirectory=argv[5];
  string RunIdentifier="BillIsAlive";
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
    slope=atof(argv[9]);
  double RMAX=76.5;
  if(argc >= 11)
    RMAX=atof(argv[10]);
  double x0=16.5;
  if(argc >= 12)
    x0=atof(argv[11]);
  double y0=17.5;
  if(argc >= 13)
    y0=atof(argv[12]);
  double z0=15.5;
  if(argc >= 14)
    z0=atof(argv[13]);
  double mm=15.5;
  if(argc >= 15)
    mm=atof(argv[14]);
  double BOXLENGTH=2.2*RMAX;
  if(argc >= 16)
    BOXLENGTH=atof(argv[15]);
  double SHRINK=1.0;
  if(argc >= 17)
    SHRINK=atof(argv[16]);
  int RANDOMSEED=8765;
  if(argc >=18)
    RANDOMSEED=atoi(argv[17]);
  int withparts=1;
  if(argc >=19)
    withparts=atoi(argv[18]);
  if(RANK == 0)
    {
      int ar=0;
      while(ar < argc)
	cerr << " " << argv[ar++] << "\n";
    }
  double BOXLENGTH2=BOXLENGTH/2.0;
  vector<double>xmin{x0-BOXLENGTH2,y0-BOXLENGTH2,z0-BOXLENGTH2};
  vector<double>xmax{x0+BOXLENGTH2,y0+BOXLENGTH2,z0+BOXLENGTH2};
  vector<double>xmini=xmin;
  vector<double>xmaxy=xmax;
  
  FractalSpace::Fractal_Memory* PFM=
    FractalSpace::FractalGravityIsolatedFirstTime(
						  TalkToMe,
						  GridLength,
						  FractalNodes0,
						  FractalNodes1,
						  FractalNodes2,
						  BaseDirectory,
						  RunIdentifier);


  FractalSpace::fractal_memory_setup(PFM);

  PFM->setNumberParticles(NumberParticles);
  PFM->level_max=10;

  for(int round=0;round < 3;round++)
    {
      FractalSpace::fractal_create(PFM);

      srand(RANDOMSEED+RANK);
      double slope3=slope+3.0;
      double expo=1.0/slope3;
      double rand_max=(double)RAND_MAX;
      double twopi=8.0*atan(1.0);
      vector<double>posx;
      vector<double>posy;
      vector<double>posz;
      vector<double>masses;
      double m0=mm/(double)(FRN*NumberParticles);

      FILE* PFPos=PFM->p_file->PFPos;

      for(int ni=0;ni<NumberParticles;ni++)
	{
	  double r0=FractalSpace::Fractal::my_rand_not_zero(rand_max);
	  double r1=RMAX*pow(r0,expo);
	  double zr=2.0*FractalSpace::Fractal::my_rand(rand_max)-1.0;
	  double Rr=sqrt(max(1.0-zr*zr,0.0));
	  double phi=twopi*FractalSpace::Fractal::my_rand(rand_max);
	  double xr=Rr*cos(phi);
	  double yr=Rr*sin(phi);
	  posx.push_back(r1*xr+x0);
	  posy.push_back(r1*yr+y0);
	  posz.push_back(r1*zr+z0);
	  masses.push_back(m0);
	}

      FractalSpace::add_particles(PFM,0,NumberParticles,posx,posy,posz,masses);
  
      FractalSpace::FractalCube(PFM,SHRINK,xmin,xmax,xmini,xmaxy);

      FractalSpace::balance_by_particles(PFM,withparts);

      FractalSpace::DoFractalGravity(PFM);

      vector<double>pot(NumberParticles);
      vector<double>accx(NumberParticles);
      vector<double>accy(NumberParticles);
      vector<double>accz(NumberParticles);

      FractalSpace::get_field(PFM,0,NumberParticles,G,xmini,xmaxy,
			      pot,accx,accy,accz);

      for(int ni=0;ni<NumberParticles;ni++)
	{
	  int lev=PFM->p_fractal->particle_list[ni]->get_highest_level();
	  double dx=posx[ni]-x0;
	  double dy=posy[ni]-y0;
	  double dz=posz[ni]-z0;
	  double r2=dx*dx+dy*dy+dz*dz+1.0e-30;
	  double r1=sqrt(r2);
	  double r3=r1*r2;
	  double MR=mm*pow(r1/RMAX,slope3);
	  double acc0=G*MR/r3;
	  double accx0=-acc0*dx;
	  double accy0=-acc0*dy;
	  double accz0=-acc0*dz;
	  acc0*=r1;
	  double acc=sqrt(pow(accx[ni],2)+pow(accy[ni],2)+pow(accz[ni],2));
	  fprintf(PFPos," %6d L %2d %13.4E %13.4E %13.4E %13.4E %13.4E %13.4E %13.4E %13.4E %13.4E %13.4E %13.4E %13.4E \n",
		  ni,lev,r1,acc,acc0,posx[ni],posy[ni],posz[ni],accx[ni],accx0,accy[ni],accy0,accz[ni],accz0);
	}
  
      FractalSpace::fractal_delete(PFM);  
    }  
  FractalSpace::fractal_memory_content_delete(PFM);

  FractalSpace::fractal_memory_delete(PFM);

  PFM=0;
  
  return 0;
}
