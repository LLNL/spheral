//---------------------------------Spheral++----------------------------------//
// PyFileIO -- A python friendly version of the FileIO interface, for use
// creating python FileIO objects.
//
// This class overrides the FileIO::read methods, since BPL has problems 
// handing back references to int, bool, and double.
//
// Created by JMO, Tue Dec 27 22:12:00 PST 2005
//----------------------------------------------------------------------------//

#include "PyFileIO.hh"
#include "Utilities/DBC.hh"

#include <string>

namespace Spheral {

//------------------------------------------------------------------------------
// Constructors.
//------------------------------------------------------------------------------
PyFileIO::
PyFileIO():
  FileIO() {
}

PyFileIO::
PyFileIO(const std::string filename, AccessType access):
  FileIO(filename, access) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
PyFileIO::
~PyFileIO() {
}

}
