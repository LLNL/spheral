//---------------------------------Spheral++----------------------------------//
// Restart
// A class templated on the object type, which requires that the object 
// provide the following methods:
//   dumpState
//   restoreState
//   label
//
// Created by JMO, Wed May 27 15:39:48 PDT 2009
//----------------------------------------------------------------------------//
#ifndef __Spheral_Restart__
#define __Spheral_Restart__

#include <string>
#include "RestartHandle.hh"

namespace Spheral {
namespace DataOutput {

template<typename Object>
class Restart: public RestartHandle {
public:
  //--------------------------- Public Interface ---------------------------//
  // Constructors.
  Restart(Object& object);

  // Destructor.
  virtual ~Restart();

  //******************************************************************************
  // Methods all restartable objects must provide.
  //******************************************************************************
  // Dump the objects state to the given file.
  virtual void dumpState(FileIOSpace::FileIO& file, const std::string& pathName) const;

  // Restore state from the given file.
  virtual void restoreState(const FileIOSpace::FileIO& file, const std::string& pathName) const;

  // Provide a label for the object type to use when writing to the file.
  virtual std::string label() const;
  //******************************************************************************

private:
  //--------------------------- Private Interface ---------------------------//
  Object* mObjectPtr;
};

}
}

#ifndef __GCCXML__
#include "RestartInline.hh"
#endif

#else

// Forward declaration.
namespace Spheral{
  namespace DataOutput {
    template<typename Object> class Restart;
  }
}

#endif
