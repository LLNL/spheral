#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
int main(int argc, char* argv[])
{
  using namespace FractalSpace;
  MPI_Init(NULL,NULL);
  bool cosmo=false;
  if(argc >= 2)
    cosmo=atoi(argv[1]) != 0;
  if(cosmo)
    fractal_nina_cosmo(argc,argv);
  else
    fractal_nina_galaxy(argc,argv);
  return 0;
}
