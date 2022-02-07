//---------------------------------Spheral++----------------------------------//
// DampedLinearSpring -- contact model based on damped linear springs
//                       Cundall & Strack Geotechnique, vol. 29, no. 1,
//                       pp. 47-65, 1979.
//----------------------------------------------------------------------------//
#ifndef __Spheral_DampedLinearSpring_hh__
#define __Spheral_DampedLinearSpring_hh__

#include "ContactModelBase.hh"
#include <string>

namespace Spheral {

template<typename Dimension> class DataBase;
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
class FileIO;

template<typename Dimension>
class DampedLinearSpring : public ContactModelBase<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;

  // Constructors.
  DampedLinearSpring(const DataBase<Dimension>& dataBase,
                     const Scalar normalSpringConstant,
                     const Scalar restitutionCoefficient);

  // Destructor.
  ~DampedLinearSpring();

  //***************************************************************************
  // Required methods from contact model
  virtual Scalar timeStep(const DataBase<Dimension>& dataBase,
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs,
                                typename Dimension::Scalar time) const override;

  virtual void   evaluateDerivatives(const Scalar time,
                                     const Scalar dt,
                                     const DataBase<Dimension>& dataBase,
                                     const State<Dimension>& state,
                                           StateDerivatives<Dimension>& derivs) const override;
  

  //****************************************************************************
  Scalar normalSpringConstant() const;
  void   normalSpringConstant(Scalar x);

  Scalar restitutionCoefficient() const;
  void   restitutionCoefficient(Scalar x);

  Scalar beta() const;
  void   beta(Scalar x);

  Scalar timeStep() const;
  void   timeStep(Scalar x);

private:
  //--------------------------- Private Interface ---------------------------//
  // No default constructor, copying, or assignment.
  DampedLinearSpring();
  DampedLinearSpring(const DampedLinearSpring&);
  DampedLinearSpring& operator=(const DampedLinearSpring&);

  Scalar mNormalSpringConstant;
  Scalar mRestitutionCoefficient;
  Scalar mBeta;
  Scalar mTimeStep;
};

}

#include "DampedLinearSpringInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class DampedLinearSpring;
}

#endif
