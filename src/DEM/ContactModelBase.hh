//---------------------------------Spheral++----------------------------------//
// ContactModelBase -- base for all DEM contact Models package for Spheral++.
//----------------------------------------------------------------------------//
#ifndef __Spheral_ContactModelBase_hh__
#define __Spheral_ContactModelBase_hh__

#include <string>

namespace Spheral {

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
//template<typename Dimension, typename DataType> class Field;
//template<typename Dimension, typename DataType> class FieldList;
class FileIO;

template<typename Dimension>
class ContactModelBase {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;

  // Constructors.
  ContactModelBase();

  // Destructor.
  virtual ~ContactModelBase();

  //***************************************************************************
  // Required methods from contact model

  virtual Scalar timeStep( ) const = 0;

  virtual Vector force(const Scalar mi, const Scalar mj,
                       const Vector ri, const Vector rj,
                       const Vector vi, const Vector vj,
                       const Scalar hi, const Scalar hj) const;
  
  virtual Vector torque(const Scalar mi, const Scalar mj,
                        const Vector ri, const Vector rj,
                        const Vector vi, const Vector vj,
                        const Scalar hi, const Scalar hj) const;

};

}

#include "ContactModelBaseInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class ContactModelBase;
}

#endif
