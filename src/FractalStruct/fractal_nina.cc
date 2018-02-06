#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
int main(int argc, char* argv[])
{
  using namespace FractalSpace;
  MPI_Init(NULL,NULL);
  if(argc >= 2 && atoi(argv[1]) != 0);
    fractal_nina_cosmo(argc,argv);
  else
    fractal_nina_galaxy(argc,argv);
  return 0;
}
