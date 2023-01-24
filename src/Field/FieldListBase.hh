//---------------------------------Spheral++----------------------------------//
// FieldListBase -- A base class to provide a generic handle on FieldLists.
//
// Created by JMO, Sun Apr 11 14:28:59 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral_FieldListBase_hh__
#define __Spheral_FieldListBase_hh__

#include <vector>

namespace Spheral {

template<typename Dimension> class FieldBase;

template<typename Dimension>
class FieldListBase {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef FieldBase<Dimension>* ElementType;
  typedef FieldBase<Dimension>* value_type;    // STL compatibility
  typedef std::vector<ElementType> StorageType;

  // Constructors.
  FieldListBase();
  FieldListBase(const FieldListBase& rhs);

  // Destructor.
  virtual ~FieldListBase();

  // Assignment operator.
  FieldListBase& operator=(const FieldListBase& rhs);

protected:
  //--------------------------- Protected Interface ---------------------------//
  mutable bool mNewCoarseNodes;
  mutable bool mNewRefineNodes;

  // Provide methods for the FieldList to register with its member Fields.
  void registerWithField(const FieldBase<Dimension>& fieldBase) const;
  void unregisterFromField(const FieldBase<Dimension>& fieldBase) const;

private:
  //--------------------------- Private Interface ---------------------------//
};

}

#include "FieldListBaseInline.hh"

#else

// Forward declare the FieldListBase class.
namespace Spheral {
  class FieldListBase;
}

#endif
