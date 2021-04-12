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
  //typedef typename Dimension::Tensor Tensor;
  //typedef typename Dimension::SymTensor SymTensor;

  //typedef typename Physics<Dimension>::TimeStepType TimeStepType;

  // Constructors.
  ContactModelBase();

  // Destructor.
  virtual ~ContactModelBase();

  //***************************************************************************
  // Required methods from contact model
  virtual Scalar timeStep(const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs) const = 0;

  virtual void force(const State<Dimension>& state,
                     StateDerivatives<Dimension>& derivs) const;

  virtual void torque(const State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs) const;
  //****************************************************************************

//protected:
  // The restart registration.
  //RestartRegistrationType mRestart;

//private:
  //--------------------------- Private Interface ---------------------------//
  // No default constructor, copying, or assignment.
  //ContactModelBase();
  //ContactModelBase(const ContactModelBase&);
  //ContactModelBase& operator=(const ContactModelBase&);
};

}

#include "ContactModelBaseInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class ContactModelBase;
}

#endif
