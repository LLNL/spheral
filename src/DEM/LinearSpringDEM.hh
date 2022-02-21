//---------------------------------Spheral++----------------------------------//
// LinearSpringDEM -- contact model based on damped linear springs
//                       Cundall & Strack Geotechnique, vol. 29, no. 1,
//                       pp. 47-65, 1979.
//----------------------------------------------------------------------------//
#ifndef __Spheral_LinearSpringDEM_hh__
#define __Spheral_LinearSpringDEM_hh__

#include "Geometry/Dimension.hh"
#include "DEM/DEMDimension.hh"
#include "DEMBase.hh"
#include <string>

namespace Spheral {

template<typename Dimension> class DataBase;
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
class FileIO;

template<typename Dimension>
class LinearSpringDEM : public DEMBase<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;

  typedef typename DEMBase<Dimension>::TimeStepType TimeStepType;

  LinearSpringDEM(const DataBase<Dimension>& dataBase,
                  const Scalar normalSpringConstant,
                  const Scalar restitutionCoefficient,
                  const Scalar stepsPerCollision,
                  const Vector& xmin,
                  const Vector& xmax);

  ~LinearSpringDEM();
 
  virtual TimeStepType dt(const DataBase<Dimension>& dataBase,
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs,
                          const Scalar time) const override;

  virtual void evaluateDerivatives(const Scalar time,
                                   const Scalar dt,
                                   const DataBase<Dimension>& dataBase,
                                   const State<Dimension>& state,
                                         StateDerivatives<Dimension>& derivs) const override;
  
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
  Scalar mNormalSpringConstant;
  Scalar mRestitutionCoefficient;
  Scalar mBeta;
  Scalar mTimeStep;

  // No default constructor, copying, or assignment.
  LinearSpringDEM();
  LinearSpringDEM(const LinearSpringDEM&);
  LinearSpringDEM& operator=(const LinearSpringDEM&);
};

}

#include "LinearSpringDEMInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class LinearSpringDEM;
}

#endif
