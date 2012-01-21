//---------------------------------Spheral++----------------------------------//
// FieldListBase -- A base class to provide a generic handle on FieldLists.
//
// Created by JMO, Sun Apr 11 14:28:59 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral_FieldListBase_hh__
#define __Spheral_FieldListBase_hh__

namespace Spheral {
namespace FieldSpace {

template<typename Dimension> class FieldBase;

class FieldListBase {

public:
  //--------------------------- Public Interface ---------------------------//
  enum FieldStorageType {
    Reference = 0,
    Copy = 1
  };

  // Constructors.
  FieldListBase();
  FieldListBase(const FieldListBase& rhs);

  // Destructor.
  virtual ~FieldListBase();

  // Assignment operator.
  FieldListBase& operator=(const FieldListBase& rhs);

protected:
  //--------------------------- Private Interface ---------------------------//
  mutable bool mNewCoarseNodes;
  mutable bool mNewRefineNodes;

  // Provide methods for the FieldList to register with its member Fields.
  template<typename Dimension> void registerWithField(const FieldBase<Dimension>& fieldBase) const;
  template<typename Dimension> void unregisterFromField(const FieldBase<Dimension>& fieldBase) const;

private:
  //--------------------------- Private Interface ---------------------------//
};

}
}

#ifndef __GCCXML__
#include "FieldListBaseInline.hh"
#endif

#else

// Forward declare the FieldListBase class.
namespace Spheral {
  namespace FieldSpace {
    class FieldListBase;
  }
}

#endif
