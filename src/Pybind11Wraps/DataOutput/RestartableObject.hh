//---------------------------------Spheral++----------------------------------//
// RestartableObject
//
// This is an object that handles registering with the RestartRegistrar.
// RestartableObject is intended as a helper for python classes which want to
// be restartable.  They can simply do the following in their __init__ method:
//   self.restart = RestartableObject(self)
// and they must provide the "label", "dumpState", and "restoreState" methods
// as usual.
// This is intended solely for use making Python objects play in our restart
// setup -- please use the C++ centric methods in registerWithRestart.hh
// for C++ restarting.
//
// Created by JMO, Thu May 28 17:47:48 PDT 2009
//----------------------------------------------------------------------------//
#ifndef __Spheral_RestartableObject__
#define __Spheral_RestartableObject__

#include <string>

#include "DataOutput/registerWithRestart.hh"

// Forward declarations.
namespace Spheral {
  namespace FileIOSpace {
    class FileIO;
  }
}

namespace Spheral {
namespace DataOutput {

class RestartableObject {

public:
  //------------------------===== Public Interface =====-----------------------//
  RestartableObject(const unsigned priority);
  virtual ~RestartableObject();

  // The methods we are providing for restart.
  virtual std::string label() const = 0;
  virtual void dumpState(FileIOSpace::FileIO& file, const std::string pathName) const = 0;
  virtual void restoreState(const FileIOSpace::FileIO& file, const std::string pathName) = 0;

private:
  //-----------------------===== Private Interface =====-----------------------//
  // The restart registration.
  RestartRegistrationType mRestart;
};

}
}

#endif

