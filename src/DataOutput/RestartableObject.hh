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
#ifndef __Spheral_RestarbableObject__
#define __Spheral_RestarbableObject__

#include <string>

#ifndef __GCCXML__
#include "DataOutput/registerWithRestart.hh"
#else
#include "fakestl.hh"
#endif

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
  RestartableObject();
  RestartableObject(const unsigned priority);
  virtual ~RestartableObject();

  // All descendents *must* provide the following.
  virtual std::string label() const = 0;
  virtual void dumpState(FileIOSpace::FileIO& file, const std::string pathName) const = 0;
  virtual void restoreState(const FileIOSpace::FileIO& file, const std::string pathName) = 0;

private:
  //-----------------------===== Private Interface =====-----------------------//
#ifndef __GCCXML__
  // The restart registration.
  DataOutput::RestartRegistrationType mRestart;
#endif
};

}
}

#endif

