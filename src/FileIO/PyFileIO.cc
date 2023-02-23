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

//------------------------------------------------------------------------------
// write vector<char> as py::bytes
//------------------------------------------------------------------------------
void
PyFileIO::
write_vector_char(const std::vector<char>& value, const std::string path) {
  py::bytes buf(&(*value.begin()), value.size());
  this->write_object(buf, path);
}

//------------------------------------------------------------------------------
// read vector<char>
//------------------------------------------------------------------------------
std::vector<char>
PyFileIO::
read_vector_char(const std::string path) const {
  py::bytes buf = this->read_object(path);
  std::string bufstr = std::string(buf);
  vector<char> result(bufstr.begin(), bufstr.end());
  ENSURE(result.size() == py::len(buf));
  return result;
}

}
