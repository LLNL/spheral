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

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
//template<typename Dimension, typename DataType> class Field;
//template<typename Dimension, typename DataType> class FieldList;
class FileIO;

template<typename Dimension>
class DampedLinearSpring : public ContactModelBase<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  //typedef typename Dimension::Tensor Tensor;
  //typedef typename Dimension::SymTensor SymTensor;

  // Constructors.
  DampedLinearSpring(Scalar YoungsModulus,
                     Scalar restitutionCoefficient);

  // Destructor.
  ~DampedLinearSpring();

  //***************************************************************************
  // Required methods from contact model

  virtual Scalar timeStep() const  override;

  virtual Vector force(const Scalar mi, const Scalar mj,
                       const Vector ri, const Vector rj,
                       const Vector vi, const Vector vj,
                       const Scalar hi, const Scalar hj) const override;
  
  virtual Vector torque(const Scalar mi, const Scalar mj,
                        const Vector ri, const Vector rj,
                        const Vector vi, const Vector vj,
                        const Scalar hi, const Scalar hj) const override;

  //****************************************************************************

//protected:
  // The restart registration.
  //RestartRegistrationType mRestart;

private:
  //--------------------------- Private Interface ---------------------------//
  // No default constructor, copying, or assignment.
  DampedLinearSpring();
  DampedLinearSpring(const DampedLinearSpring&);
  DampedLinearSpring& operator=(const DampedLinearSpring&);

  Scalar mYoungsModulus;
  Scalar mRestitutionCoefficient;
  Scalar mBeta;
};

}

#include "DampedLinearSpringInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class DampedLinearSpring;
}

#endif
