//---------------------------------Spheral++----------------------------------//
// JohnsonCookStrength -- Implements the Johnson-Cook strength model.
//   ** Need reference **
//
// Created by JMO, Tue Jul 28 13:39:50 PDT 2015
//----------------------------------------------------------------------------//
#ifndef __Spheral_JohnsonCookStrength_hh__
#define __Spheral_JohnsonCookStrength_hh__

#include "StrengthModel.hh"
#include "PolynomialFit.hh"

namespace Spheral {
  namespace SolidMaterial {
    template<typename Dimension> class SolidEquationOfState;
  }
}

namespace Spheral {
namespace SolidMaterial {

template<typename Dimension>
class JohnsonCookStrength: public StrengthModel<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;

  // Constructors, destructor.
  JohnsonCookStrength(const SolidEquationOfState<Dimension>& eos,
                      const StrengthModel<Dimension>& shearModulusModel,
                      const double A,
                      const double B,
                      const double C,
                      const double C4,
                      const double m,
                      const double nhard,
                      const double epsdot0,
                      const double epsdotmin,
                      const double Tmelt,
                      const double Troom);
  virtual ~JohnsonCookStrength();

  // Override the required generic interface.
  virtual void shearModulus(FieldSpace::Field<Dimension, Scalar>& shearModulus,
                            const FieldSpace::Field<Dimension, Scalar>& density,
                            const FieldSpace::Field<Dimension, Scalar>& specificThermalEnergy,
                            const FieldSpace::Field<Dimension, Scalar>& pressure) const;

  virtual void yieldStrength(FieldSpace::Field<Dimension, Scalar>& yieldStrength,
                             const FieldSpace::Field<Dimension, Scalar>& density,
                             const FieldSpace::Field<Dimension, Scalar>& specificThermalEnergy,
                             const FieldSpace::Field<Dimension, Scalar>& pressure,
                             const FieldSpace::Field<Dimension, Scalar>& plasticStrain,
                             const FieldSpace::Field<Dimension, Scalar>& plasticStrainRate) const;

  virtual void soundSpeed(FieldSpace::Field<Dimension, Scalar>& soundSpeed,
                          const FieldSpace::Field<Dimension, Scalar>& density,
                          const FieldSpace::Field<Dimension, Scalar>& specificThermalEnergy,
                          const FieldSpace::Field<Dimension, Scalar>& pressure,
                          const FieldSpace::Field<Dimension, Scalar>& fluidSoundSpeed) const;

  // Access the strength parameters.
  double A() const;
  double B() const;
  double C() const;
  double C4() const;
  double m() const;
  double nhard() const;
  double epsdot0() const;
  double epsdotmin() const;
  double Tmelt() const;
  double Troom() const;

private:
  //--------------------------- Private Interface ---------------------------//
  const SolidEquationOfState<Dimension>* mEOSPtr;
  const StrengthModel<Dimension>* mShearModulusModelPtr;
  const double mA, mB, mC, mC4, mm, mnhard, mEpsdot0, mEpsdotmin, mTmelt, mTroom;

  // No copying or assignment.
  JohnsonCookStrength(const JohnsonCookStrength&);
  JohnsonCookStrength& operator=(const JohnsonCookStrength&);
};

}
}

#else

// Forward declaration.
namespace Spheral {
  namespace SolidMaterial {
    template<typename Dimension> class JohnsonCookStrength;
  }
}

#endif

