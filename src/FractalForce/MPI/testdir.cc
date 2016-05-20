#include <cstdlib>
#include <iostream>
#include <fstream>
#include "fftw3-mpi.h"
int main()
{
  using namespace std;
  MPI_Init(NULL, NULL);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  fftw_mpi_init();
  fftw_mpi_cleanup();

  MPI_Finalize();

  printf("Rank %d done\n", rank);  fflush(stdout);

  return 0;
}
