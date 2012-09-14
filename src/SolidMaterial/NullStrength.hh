//---------------------------------Spheral++----------------------------------//
// NullStrength -- An implentation of StrengthModel which mimics having
// no strength.
//
// Created by JMO, Thu Sep 16 23:48:47 PDT 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral_NullStrength_hh__
#define __Spheral_NullStrength_hh__

#include "StrengthModel.hh"

namespace Spheral {
namespace SolidMaterial {

template<typename Dimension>
class NullStrength: public StrengthModel<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Constructors, destructor.
  NullStrength() {};
  virtual ~NullStrength() {};

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

  // No copying or assignment.
  NullStrength(const NullStrength&);
  NullStrength& operator=(const NullStrength&);
};

#ifndef __GCCXML__

//------------------------------------------------------------------------------
// Compute the shear modulus.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
NullStrength<Dimension>::
shearModulus(const double density,
             const double specificThermalEnergy,
             const double pressure) const {
  return 0.0;
}

//------------------------------------------------------------------------------
// Compute the yeild strength.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
NullStrength<Dimension>::
yieldStrength(const double density,
              const double specificThermalEnergy,
              const double pressure,
              const double plasticStrain,
              const double plasticStrainRate) const {
  return 0.0;
}

//------------------------------------------------------------------------------
// Compute the full sound speed.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
NullStrength<Dimension>::
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
    template<typename Dimension> class NullStrength;
  }
}

#endif

