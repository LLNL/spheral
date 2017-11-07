#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
int main()
{
  using namespace FractalSpace;
  cout << "starting out " << endl;
  Fractal_Memory* p_fractal_memory= new Fractal_Memory;
  //  fractal_memory_parameters(*p_fractal_memory);  
  bool MR=p_fractal_memory->MPIrun;
  int GR=p_fractal_memory->grid_length;
  bool PR=p_fractal_memory->periodic;
  int NP=p_fractal_memory->number_particles;
  Mess* p_mess=new Mess(MR,GR,PR,NP);
  p_fractal_memory->p_mess=p_mess;
  string BD=p_fractal_memory->BaseDirectory;
  int FR=p_mess->FractalRank;
  string RUN=p_fractal_memory->RUN;
  File* p_file=new File(BD,FR,RUN);
  p_fractal_memory->p_file=p_file;
  //  fractal_force_init(p_fractal_memory);
  Fractal* p_fractal=new Fractal(*p_fractal_memory);
  //  p_fractal->p_mess=p_mess;
  //  p_fractal->p_file=p_file;
  cout << " start " << p_fractal_memory << " " << p_fractal << endl;
  fractal_force(*p_fractal,*p_fractal_memory);
  delete p_mess;
  p_mess=0;
  delete p_fractal_memory;
  p_fractal_memory=0;
  delete p_file;
  p_file=0;
  delete p_fractal;
  p_fractal=0;
  int result=1;
  return result;
}
