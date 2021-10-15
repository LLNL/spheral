// managed_parallel.cpp
#include <iostream>

#include "../Timer.hh"

#ifdef USE_MPI
#include <mpi.h>
#include "Distributed/Communicator.hh"
#endif

Timer Everything("Everything");
Timer loop1(" Loop 1", Everything);
Timer loop2(" Loop 2", Everything);

int main(int argc, char **argv) {

  int rank, number_procs;
#ifdef USE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(Spheral::Communicator::communicator(), &rank);
  MPI_Comm_size(Spheral::Communicator::communicator(), &number_procs);
#else
  rank=0;
  number_procs=1;
#endif
  cout << "rank = " << rank << " number_procs= " << number_procs << endl;


  Everything.start();
  double y1=0, y2=0;

  loop1.start();
  for(int i=0; i<(rank+1)*12345678; i++) {
    y1 += (float)i;
  }
  loop1.stop();

  loop2.start();
  for(int i=0; i<(rank+1)*12345678; i++) {
    y2 += (float)i/(i+1);
  }
  loop2.stop();

  cout << " y1= " << y1 << " y2= " << y2 << endl;
  Everything.stop();

  Timer::TimerSummary(rank, number_procs);

#ifdef USE_MPI
  MPI_Finalize();
#endif
}

