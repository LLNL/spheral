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
#else
typedef int MPI_Comm;
#endif

namespace Spheral {

class Communicator {

public:
  //------------------------===== Public Interface =====-----------------------//
  // Get the instance.
  static Communicator& instance() { static Communicator theInstance; return theInstance; }

  // Access the communicator.
  static MPI_Comm& communicator() { return instance().mCommunicator; }
  static void communicator(MPI_Comm& comm) { instance().mCommunicator = comm; }
  static MPI_Comm* comm_ptr() {
#ifdef USE_MPI
    return &(instance().mCommunicator);
#else
    return nullptr;
#endif
  }
  static void finalize() {
#ifdef USE_MPI
    int finalize = MPI_Finalize();
    VERIFY(finalize);
#endif
  }

private:
  //------------------------===== Private Interface =====----------------------//
#ifdef USE_MPI
  MPI_Comm mCommunicator;
#else
  MPI_Comm mCommunicator = 0;
#endif

  // No public constructors, destructor, or assignment.
  Communicator();
  Communicator(const Communicator&);
  Communicator& operator=(const Communicator&);
  ~Communicator();
};

}

#endif
