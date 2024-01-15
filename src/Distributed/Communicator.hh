//---------------------------------Spheral++----------------------------------//
// Communicator
// A singleton object which holds the MPI communicator Spheral uses.
//
// Created by JMO, Wed Jan 30 14:52:05 PST 2013
//----------------------------------------------------------------------------//
#ifndef __Spheral_Communicator__
#define __Spheral_Communicator__

#ifdef USE_MPI
#include <mpi.h>
#endif

namespace Spheral {

class Communicator {

public:
  //------------------------===== Public Interface =====-----------------------//
  // Get the instance.
  static Communicator& instance() { static Communicator theInstance; return theInstance; }

  // Access the communicator.
#ifdef USE_MPI
  static MPI_Comm& communicator() { return instance().mCommunicator; }
  static void communicator(MPI_Comm& comm) { instance().mCommunicator = comm; }
#else
  static int communicator() { return 0; }
  static void communicator(int&) {}
#endif

private:
  //------------------------===== Private Interface =====----------------------//
#ifdef USE_MPI
  MPI_Comm mCommunicator;
#endif

  // No public constructors, destructor, or assignment.
  Communicator();
  Communicator(const Communicator&);
  Communicator& operator=(const Communicator&);
  ~Communicator();
};

}

#endif
