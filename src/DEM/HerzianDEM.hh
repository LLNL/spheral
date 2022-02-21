//---------------------------------Spheral++----------------------------------//
// HerzianDEM -- contact model based on damped linear springs
//                       Cundall & Strack Geotechnique, vol. 29, no. 1,
//                       pp. 47-65, 1979.
//----------------------------------------------------------------------------//
#ifndef __Spheral_HerzianDEM_hh__
#define __Spheral_HerzianDEM_hh__

#include "DEMBase.hh"

namespace Spheral {

template<typename Dimension> class DataBase;
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
class FileIO;

template<typename Dimension>
class HerzianDEM : public DEMBase<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;

  typedef typename DEMBase<Dimension>::TimeStepType TimeStepType;
  
  HerzianDEM(const DataBase<Dimension>& dataBase,
             const Scalar YoungsModulus,
             const Scalar restitutionCoefficient,
             const Scalar stepsPerCollision,
             const Vector& xmin,
             const Vector& xmax);

  ~HerzianDEM();

  virtual TimeStepType dt(const DataBase<Dimension>& dataBase,
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs,
                          const Scalar time) const override;

  virtual void   evaluateDerivatives(const Scalar time,
                                     const Scalar dt,
                                     const DataBase<Dimension>& dataBase,
                                     const State<Dimension>& state,
                                           StateDerivatives<Dimension>& derivs) const override;
  

  //****************************************************************************
  Scalar YoungsModulus() const;
  void   YoungsModulus(Scalar x);

  Scalar restitutionCoefficient() const;
  void   restitutionCoefficient(Scalar x);

  Scalar beta() const;
  void   beta(Scalar x);

private:
  Scalar mYoungsModulus;
  Scalar mRestitutionCoefficient;
  Scalar mBeta;

  //--------------------------- Private Interface ---------------------------//
  // No default constructor, copying, or assignment.
  HerzianDEM();
  HerzianDEM(const HerzianDEM&);
  HerzianDEM& operator=(const HerzianDEM&);
};

}

#include "HerzianDEMInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class HerzianDEM;
}

#endif
