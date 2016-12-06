#include <cstdlib>
#include <iostream>
#include <fstream>
#include "fftw3-mpi.h"
int main()
{
  using namespace std;
  MPI::Init();
  fftw_mpi_init();
  fftw_mpi_cleanup();
  MPI::Finalize();
  return 0;
}
