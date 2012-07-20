//---------------------------------Spheral++----------------------------------//
// StrengthModel -- The interface base class for strength models.
//
// Created by JMO, Wed Sep 8 15:18:44 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral_StrengthModel_hh__
#define __Spheral_StrengthModel_hh__

#include "Infrastructure/SpheralFunctions.hh"
#include "Utilities/DBC.hh"

namespace Spheral {
namespace SolidMaterial {

class StrengthModel {
public:
  //--------------------------- Public Interface ---------------------------//
  // Constructors, destructor.
  StrengthModel() {};
  virtual ~StrengthModel() {};

  // The generic interface we require all strength models to provide.
  virtual double shearModulus(const double density,
                              const double specificThermalEnergy,
                              const double pressure) const = 0;

  virtual double yieldStrength(const double density,
                               const double specificThermalEnergy,
                               const double pressure,
                               const double plasticStrain,
                               const double plasticStrainRate) const = 0;

  // An overridable method to compute the full sound speed.
  virtual double soundSpeed(const double density,
                            const double specificThermalEnergy,
                            const double pressure,
                            const double fluidSoundSpeed) const;

private:
  //--------------------------- Private Interface ---------------------------//

  // No copying or assignment.
  StrengthModel(const StrengthModel&);
  StrengthModel& operator=(const StrengthModel&);
};

//------------------------------------------------------------------------------
// Compute the full sound speed.
//------------------------------------------------------------------------------
inline
double
StrengthModel::
soundSpeed(const double density,
           const double specificThermalEnergy,
           const double pressure,
           const double fluidSoundSpeed) const {
  REQUIRE(distinctlyGreaterThan(density, 0.0));
  const double cs2 = fluidSoundSpeed*fluidSoundSpeed + 4.0/3.0 * shearModulus(density, specificThermalEnergy, pressure) / density;
  ENSURE(cs2 > 0.0);
  return std::sqrt(cs2);
}

}
}

#else

// Forward declaration.
namespace Spheral {
  namespace SolidMaterial {
    class StrengthModel;
  }
}

#endif

