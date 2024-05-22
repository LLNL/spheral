//---------------------------------Spheral++----------------------------------//
// Communicator
// A singleton object which holds the MPI communicator Spheral uses.
//
// Created by JMO, Wed Jan 30 14:52:05 PST 2013
//----------------------------------------------------------------------------//
#include "Communicator.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Default constructor (private).
//------------------------------------------------------------------------------
Communicator::
Communicator()
#ifdef USE_MPI
  : mCommunicator() {
  mCommunicator = MPI_COMM_WORLD;
#else
{
#endif
}

//------------------------------------------------------------------------------
// Destructor (private).
//------------------------------------------------------------------------------
Communicator::
~Communicator() {
}

}
