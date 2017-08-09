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
// Get the instance.
//------------------------------------------------------------------------------
Communicator&
Communicator::
instance() {
  return *Communicator::instancePtr();
}

//------------------------------------------------------------------------------
// Get the instance (pointer).
//------------------------------------------------------------------------------
Communicator*
Communicator::
instancePtr() {
  if (mInstancePtr == 0) mInstancePtr = new Communicator;
  CHECK(mInstancePtr != 0);
  return mInstancePtr;
}

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

//------------------------------------------------------------------------------
// Initialize the static instance pointer.
//-----------------------------------------------------------------------------
Spheral::Communicator* Spheral::Communicator::mInstancePtr = 0;

