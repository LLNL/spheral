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
  Mess::IAMROOT=Ranky == 21;
  string _disK_="d";
  if(argc >= 2)
    _disK_=argv[1];
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
  int FR0=dims[0];
  int FR1=dims[1];
  int FR2=dims[2];
  int _mulT_=4;
  if(argc >= 7)
    _mulT_=atoi(argv[6]);
  int PADDING=-1;
  if(argc >= 8)
    PADDING=atoi(argv[7]);
  int node_load=20000;
  if(argc >= 9)
    node_load=atoi(argv[8]);
  if(Ranky == 0)
    {
    cout << "starting out " << argc << " " << argv[0] << " " << FRN << " " << _disK_ << " " << GRL << " " << FR0 << " " << FR1 << " " << FR2 << " " << _mulT_ << " " << PADDING << " " << node_load << "\n";
    int ar=0;
    while(ar < argc)
      cerr << " " << argv[ar++] << "\n";
    }
  Fractal_Memory* PFM= new Fractal_Memory;
  PFM->FractalNodes0=FR0;
  PFM->FractalNodes1=FR1;
  PFM->FractalNodes2=FR2;
  PFM->grid_length=GRL;
  PFM->hypre_max_node_load=node_load;
  PFM->hypre_multiplier=-1;
  fractal_memory_parameters(PFM,_disK_,_mulT_);  
  PADDING=max(-1,min(1,PADDING));
  PFM->padding=PADDING;
  // PFM->MPIrun=true;
  int GR=PFM->grid_length;
  bool PR=PFM->periodic;
  int NP=PFM->number_particles;
  int FN=PFM->FFTNodes;
  MPI_Comm FW=Fractal_Memory::FRACTAL_UNIVERSE;
  FR0=PFM->FractalNodes0;
  FR1=PFM->FractalNodes1;
  FR2=PFM->FractalNodes2;
  Mess* p_mess=new Mess(true,GR,PR,NP,FR0,FR1,FR2,FN,FW);
  PFM->p_mess=p_mess;
  PFM->FFTNodes=PFM->p_mess->FFTNodes;
  string BD=PFM->BaseDirectory;
  int FR=p_mess->FractalRank;
  string RUN=PFM->RUN;
  File* p_file=new File(BD,FRN,FR,RUN);
  PFM->p_file=p_file;
  PFM->p_mess->p_file=p_file;
  fractal_force_init(PFM);
  Fractal* p_fractal=new Fractal(*PFM);
  int result=fractal_force_wrapper(PFM,p_fractal);
  delete p_mess;
  p_mess=0;
  delete PFM;
  PFM=0;
  delete p_file;
  p_file=0;
  delete p_fractal;
  p_fractal=0;
  return result;
}
