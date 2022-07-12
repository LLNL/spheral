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
                  const Scalar normalRestitutionCoefficient,
                  const Scalar tangentialSpringConstant,
                  const Scalar tangentialRestitutionCoefficient,
                  const Scalar dynamicFrictionCoefficient,
                  const Scalar staticFrictionCoefficient,
                  const Scalar rollingFrictionCoefficient,
                  const Scalar torsionalFrictionCoefficient,
                  const Scalar cohesiveTensileStrength,
                  const Scalar shapeFactor,
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

  Scalar normalRestitutionCoefficient() const;
  void   normalRestitutionCoefficient(Scalar x);

  Scalar tangentialSpringConstant() const;
  void   tangentialSpringConstant(Scalar x);

  Scalar tangentialRestitutionCoefficient() const;
  void   tangentialRestitutionCoefficient(Scalar x);

  Scalar dynamicFrictionCoefficient() const;
  void   dynamicFrictionCoefficient(Scalar x);

  Scalar staticFrictionCoefficient() const;
  void   staticFrictionCoefficient(Scalar x);

  Scalar rollingFrictionCoefficient() const;
  void   rollingFrictionCoefficient(Scalar x);

  Scalar torsionalFrictionCoefficient() const;
  void   torsionalFrictionCoefficient(Scalar x);

  Scalar cohesiveTensileStrength() const;
  void   cohesiveTensileStrength(Scalar x);

  Scalar shapeFactor() const;
  void   shapeFactor(Scalar x);

  Scalar normalBeta() const;
  void   normalBeta(Scalar x);

  Scalar tangentialBeta() const;
  void   tangentialBeta(Scalar x);

  Scalar timeStep() const;
  void   timeStep(Scalar x);

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const override { return "LinearSpringDEM" ; }
  //****************************************************************************

private:
  //--------------------------- Private Interface ---------------------------//
  Scalar mNormalSpringConstant;
  Scalar mNormalRestitutionCoefficient;
  Scalar mTangentialSpringConstant;
  Scalar mTangentialRestitutionCoefficient;
  Scalar mDynamicFrictionCoefficient;     // coefficient of friciton - dynamic
  Scalar mStaticFrictionCoefficient;      // coefficient of friction - static
  Scalar mRollingFrictionCoefficient;     // coefficient of friction - rolling
  Scalar mTorsionalFrictionCoefficient;   // coefficient of friction - torsional 
  Scalar mCohesiveTensileStrength;
  Scalar mShapeFactor;                    // varies between 0 and 1 to account to non spherical shapes & influences rolling/torsion spring parameters

  Scalar mNormalBeta;
  Scalar mTangentialBeta;
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
