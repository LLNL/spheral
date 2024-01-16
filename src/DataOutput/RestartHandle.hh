//---------------------------------Spheral++----------------------------------//
// RestartHandle
// An untyped interface class to Restart<Object>.  Alows us to call the 
// following restart methods without knowing the type of object we're working
// with:
//   dumpState
//   restoreState
//   label
//
// Created by JMO, Wed May 27 15:39:48 PDT 2009
//----------------------------------------------------------------------------//
#ifndef __Spheral_RestartHandle__
#define __Spheral_RestartHandle__

#include <string>

// Forward declarations.
namespace Spheral {
  class FileIO;
}

namespace Spheral {

class RestartHandle {
public:
  //--------------------------- Public Interface ---------------------------//
  // Constructors.
  RestartHandle() {};

  // Destructor.
  virtual ~RestartHandle() {};

  //******************************************************************************
  // Methods all restartable objects must provide.
  //******************************************************************************
  // Dump the objects state to the given file.
  virtual void dumpState(FileIO& file, const std::string& pathName) const = 0;

  // Restore state from the given file.
  virtual void restoreState(const FileIO& file, const std::string& pathName) const = 0;

  // Provide a label for the object type to use when writing to the file.
  virtual std::string label() const = 0;
  //******************************************************************************

private:
  //--------------------------- Private Interface ---------------------------//
};

}

#endif
