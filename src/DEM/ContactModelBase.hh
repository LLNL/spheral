//---------------------------------Spheral++----------------------------------//
// ContactModelBase -- base for all DEM contact Models package for Spheral++.
//----------------------------------------------------------------------------//
#ifndef __Spheral_ContactModelBase_hh__
#define __Spheral_ContactModelBase_hh__

#include <string>

namespace Spheral {

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class DataBase;
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

  virtual Scalar timeStep(const DataBase<Dimension>& dataBase,
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs,
                          const typename Dimension::Scalar time) const = 0;

  virtual void   evaluateDerivatives(const Scalar time,
                                     const Scalar dt,
                                     const DataBase<Dimension>& dataBase,
                                     const State<Dimension>& state,
                                           StateDerivatives<Dimension>& derivs) const = 0;

};

}

#include "ContactModelBaseInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class ContactModelBase;
}

#endif
