//---------------------------------Spheral++----------------------------------//
// RestartableObject
// This is an object that handles registering it's descendents with the 
// RestartRegistrar.  The idea here is you just inherit from this bad boy, 
// and blamo you're restartable.
// This is intended solely for use making Python objects play in our restart
// setup -- please use the C++ centric methods in registerWithRestart.hh
// for C++ restarting.
//
// Created by JMO, Thu May 28 17:47:48 PDT 2009
//----------------------------------------------------------------------------//
#include "RestartableObject.hh"
#include "FileIO/FileIO.hh"

namespace Spheral {
namespace DataOutput {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
RestartableObject::
RestartableObject(const unsigned priority):
  mRestart(registerWithRestart(*this, priority)) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
RestartableObject::
~RestartableObject() {
}

}
}
