//---------------------------------Spheral++----------------------------------//
// ConstantStrength -- An implentation of StrengthModel returning constant
// values for the shear modulus and yeild strength.
//
// Created by JMO, Mon Jul 24 10:34:40 PDT 2006
//----------------------------------------------------------------------------//
#ifndef __Spheral_ConstantStrength_hh__
#define __Spheral_ConstantStrength_hh__

#include "StrengthModel.hh"

namespace Spheral {
namespace SolidMaterial {

class ConstantStrength: public StrengthModel {
public:
  //--------------------------- Public Interface ---------------------------//
  // Constructors, destructor.
  ConstantStrength(const double mu0,
                   const double Y0);
  virtual ~ConstantStrength() {};

  // The generic interface we require all strength models to provide.
  virtual double shearModulus(const double density,
                              const double specificThermalEnergy,
                              const double pressure) const;

  virtual double yieldStrength(const double density,
                               const double specificThermalEnergy,
                               const double pressure,
                               const double plasticStrain,
                               const double plasticStrainRate) const;

  // An overridable method to compute the full sound speed.
  virtual double soundSpeed(const double density,
                            const double specificThermalEnergy,
                            const double pressure,
                            const double fluidSoundSpeed) const;

private:
  //--------------------------- Private Interface ---------------------------//
  // The values for shear modulus and yeild strength.
  double mShearModulus0;
  double mYeildStrength0;

  // No copying or assignment.
  ConstantStrength(const ConstantStrength&);
  ConstantStrength& operator=(const ConstantStrength&);
};

#ifndef __GCCXML__

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
inline
ConstantStrength::
ConstantStrength(const double mu0,
                 const double Y0):
  mShearModulus0(mu0),
  mYeildStrength0(Y0) {
}

//------------------------------------------------------------------------------
// Compute the shear modulus.
//------------------------------------------------------------------------------
inline
double
ConstantStrength::
shearModulus(const double density,
             const double specificThermalEnergy,
             const double pressure) const {
  return mShearModulus0;
}

//------------------------------------------------------------------------------
// Compute the yeild strength.
//------------------------------------------------------------------------------
inline
double
ConstantStrength::
yieldStrength(const double density,
              const double specificThermalEnergy,
              const double pressure,
              const double plasticStrain,
              const double plasticStrainRate) const {
  return mYeildStrength0;
}

//------------------------------------------------------------------------------
// Compute the full sound speed.
//------------------------------------------------------------------------------
inline
double
ConstantStrength::
soundSpeed(const double density,
           const double specificThermalEnergy,
           const double pressure,
           const double fluidSoundSpeed) const {
  return fluidSoundSpeed;
}

#endif

}
}

#else

// Forward declaration.
namespace Spheral {
  namespace SolidMaterial {
    class ConstantStrength;
  }
}

#endif

