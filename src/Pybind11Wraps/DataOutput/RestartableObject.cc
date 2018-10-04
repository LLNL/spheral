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

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
RestartableObject::
RestartableObject(py::object& self, const unsigned priority):
  mRestart(registerWithRestart(*this, priority)),
  mSelf(self) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
RestartableObject::
~RestartableObject() {
}

//------------------------------------------------------------------------------
// label
//------------------------------------------------------------------------------
std::string
RestartableObject::
label() const {
  py::object self_label = mSelf.attr("label");
  return py::str(self_label());
}

//------------------------------------------------------------------------------
// dumpState
//------------------------------------------------------------------------------
void
RestartableObject::
dumpState(FileIO& fileIO, const std::string& pathName) const {
  py::object self_dumpState = mSelf.attr("dumpState");
  self_dumpState(fileIO, pathName);
}

//------------------------------------------------------------------------------
// restoreState
//------------------------------------------------------------------------------
void
RestartableObject::
restoreState(const FileIO& fileIO, const std::string& pathName) {
  py::object self_restoreState = mSelf.attr("restoreState");
  self_restoreState(fileIO, pathName);
}

}
