#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
int main(int argc, char* argv[])
{
  using namespace FractalSpace;
  MPI_Init(NULL,NULL);
  int FRN;
  MPI_Comm_size(MPI_COMM_WORLD,&FRN);
  //
  // Intel/IBM (1/0) default Intel
  // Gridlength      default 256
  // dims[0]         default 0
  // dims[1]         default 0
  // dims[2]         default 0
  // _mulT_          default 4
  //
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
  int FR0=dims[0];
  int FR1=dims[1];
  int FR2=dims[2];
  int _mulT_=4;
  if(argc >= 7)
    _mulT_=atoi(argv[6]);
  cout << "starting out " << argc << " " << FRN << " " << _inteL_ << " " << GRL << " " << FR0 << " " << FR1 << " " << FR2;
  cout << " " << _mulT_ << "\n";
  Fractal_Memory* p_fractal_memory= new Fractal_Memory;
  p_fractal_memory->FractalNodes0=FR0;
  p_fractal_memory->FractalNodes1=FR1;
  p_fractal_memory->FractalNodes2=FR2;
  p_fractal_memory->grid_length=GRL;
  fractal_memory_parameters(*p_fractal_memory,_inteL_,_mulT_);  
  int FractalRank;
  //  MPI_Comm_rank(MPI_COMM_WORLD,&FractalRank);
  vector <int> BoxA;
  //  BoxA=p_fractal_memory->Boxes[FractalRank];
  bool MR=p_fractal_memory->MPIrun;
  int GR=p_fractal_memory->grid_length;
  bool PR=p_fractal_memory->periodic;
  int NP=p_fractal_memory->number_particles;
  int FN=p_fractal_memory->FFTNodes;
  MPI_Comm FW=MPI_COMM_WORLD;
  FR0=p_fractal_memory->FractalNodes0;
  FR1=p_fractal_memory->FractalNodes1;
  FR2=p_fractal_memory->FractalNodes2;
  Mess* p_mess=new Mess(MR,GR,PR,NP,FR0,FR1,FR2,FN,FW);
  p_fractal_memory->p_mess=p_mess;
  p_fractal_memory->FFTNodes=p_fractal_memory->p_mess->FFTNodes;
  string BD=p_fractal_memory->BaseDirectory;
  int FR=p_mess->FractalRank;
  string RUN=p_fractal_memory->RUN;
  File* p_file=new File(BD,FRN,FR,RUN);
  //  p_file->FileFractal << " before A " << BoxA[0] << " " << BoxA[1] << " " << BoxA[2] << " " << BoxA[3] << " " << BoxA[4] << " " << BoxA[5] << "\n";
  p_fractal_memory->p_file=p_file;
  p_fractal_memory->p_mess->p_file=p_file;
  FractalRank=p_mess->FractalRank;
  fractal_force_init(p_fractal_memory);
  BoxA=p_fractal_memory->Boxes[FractalRank];
  p_file->FileFractal << " before B " << BoxA[0] << " " << BoxA[1] << " " << BoxA[2] << " " << BoxA[3] << " " << BoxA[4] << " " << BoxA[5] << "\n";
  Fractal* p_fractal=new Fractal(*p_fractal_memory);
  BoxA=p_fractal_memory->Boxes[FractalRank];
  vector <int> BoxB;
  p_fractal->getBox(BoxB);
  p_file->FileFractal << " before C " << BoxA[0] << " " << BoxA[1] << " " << BoxA[2] << " " << BoxA[3] << " " << BoxA[4] << " " << BoxA[5] << "\n";
  p_file->FileFractal << " before D " << BoxB[0] << " " << BoxB[1] << " " << BoxB[2] << " " << BoxB[3] << " " << BoxB[4] << " " << BoxB[5] << "\n";
  int result=fractal_force_wrapper(p_fractal_memory,p_fractal);
  delete p_mess;
  p_mess=0;
  delete p_fractal_memory;
  p_fractal_memory=0;
  delete p_file;
  p_file=0;
  delete p_fractal;
  p_fractal=0;
  return result;
}
