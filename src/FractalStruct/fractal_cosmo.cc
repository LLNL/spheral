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
  Fractal_Memory::FRACTAL_UNIVERSE=MPI_COMM_WORLD;
  int FRN;
  MPI_Comm_size(Fractal_Memory::FRACTAL_UNIVERSE,&FRN);
  int Ranky;
  MPI_Comm_rank(Fractal_Memory::FRACTAL_UNIVERSE,&Ranky);
  Mess::IAMROOT=Ranky == 0;
  /*
    BaseDirectory
    RUN
    GridLength
    FractalNodes0
    FractalNodes1
    FractalNodes2
    Multiplier
    PADDING
    NodeLoad
  */
  string BaseDirectory="/p/lscratche/jensv/cosmo/";
  if(argc >= 2)
    BaseDirectory=argv[1];
  string RUN="BillLives";
  if(argc >= 3)
    RUN=argv[2];
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
  double _mulT_=4.0;
  if(argc >= 8)
    _mulT_=atof(argv[7]);
  int PADDING=-1;
  if(argc >= 9)
    PADDING=atoi(argv[8]);
  int node_load=20000;
  if(argc >= 10)
    node_load=atoi(argv[9]);
  if(Ranky == 0)
    {
      int ar=0;
      while(ar < argc)
	cerr << " " << argv[ar++] << "\n";
    }
  Fractal_Memory* PFM= fractal_memory_create();
  PFM->FractalNodes0=dims[0];
  PFM->FractalNodes1=dims[1];
  PFM->FractalNodes2=dims[2];
  PFM->FractalNodes=FRN;
  PFM->grid_length=GRL;
  PFM->hypre_max_node_load=node_load;
  PFM->BaseDirectory=BaseDirectory;
  PFM->RUN=RUN;
  fractal_memory_parameters(PFM,_mulT_);  
  PADDING=max(-1,min(1,PADDING));
  PFM->padding=PADDING;
  PFM->periodic=true;
  Mess* p_mess=new Mess(true,
			PFM->grid_length,
			true,
			PFM->number_particles,
			PFM->FractalNodes0,
			PFM->FractalNodes1,
			PFM->FractalNodes2,
			PFM->FFTNodes,
			Fractal_Memory::FRACTAL_UNIVERSE);
  PFM->p_mess=p_mess;
  PFM->FFTNodes=PFM->p_mess->FFTNodes;
  File* p_file=0;
  if(BaseDirectory == "")
    p_file=new File();
  else
    p_file=new File(PFM->BaseDirectory,
		    PFM->FractalNodes,
		    PFM->p_mess->FractalRank,
		    PFM->RUN);
  PFM->p_file=p_file;
  PFM->p_mess->p_file=p_file;
  fractal_force_init(PFM);
  Fractal* PF=new Fractal(*PFM);
  int result=fractal_force_wrapper(PFM,PF);
  fractal_memory_content_delete(PFM);
  fractal_memory_delete(PFM);
  PFM=0;
  return result;
}
