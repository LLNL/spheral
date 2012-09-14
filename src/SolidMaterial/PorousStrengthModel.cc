//---------------------------------Spheral++----------------------------------//
// PorousStrengthModel
//
// An implementation of strain-alpha porosity model described in
// Wunnemann, Collins, & Melosh, Icarus, 180, 514-527 (2006)
// "A strain-based porosity model for use in hydrocode simulations of impacts
//  and implications for transient crater growth in porous targets"
//
// This model assumes you will provide a solid EOS which will be modified.
// The underlying actualy solid EOS should provide the reference density, which
// will be treated here as the compacted true solid reference density.
//
// Note this model introduces a new state variable, the distention (alpha), which
// the pressure now depends on.  This implies our usual definition of P(rho, eps)
// now becomes P(rho, eps, alpha).  Our EOS interface does not recognize this
// this parameter, so we store alpha locally and only allow Field updates of the
// pressure (forbidding the single value P lookup the EOS usually allows).
//
// Created by JMO, Thu Sep 13 15:50:28 PDT 2012
//----------------------------------------------------------------------------//
#include "PorousStrengthModel.hh"

namespace Spheral {
namespace SolidMaterial {

using namespace std;
using std::abs;
using std::min;
using std::max;

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
PorousStrengthModel<Dimension>::
PorousStrengthModel(const StrengthModel& strengthModel):
  mSolidStrength(strengthModel),
  mAlphaPtr(0) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
PorousStrengthModel<Dimension>::
~PorousStrengthModel() {
}

//------------------------------------------------------------------------------
// Compute the shear modulus.
//------------------------------------------------------------------------------
template<typename Dimension>
double
PorousStrengthModel<Dimension>::
shearModulus(const double density,
             const double specificThermalEnergy,
             const double pressure) const {

  // Find the solid shear modulus.
  

}

//------------------------------------------------------------------------------
// Compute the yield strength.
//------------------------------------------------------------------------------
template<typename Dimension>
double
PorousStrengthModel<Dimension>::
yieldStrength(const double density,
              const double specificThermalEnergy,
              const double pressure,
              const double plasticStrain,
              const double plasticStrainRate) const {
}

//------------------------------------------------------------------------------
// Compute the full sound speed.
//------------------------------------------------------------------------------
template<typename Dimension>
double
PorousStrengthModel<Dimension>::
soundSpeed(const double density,
           const double specificThermalEnergy,
           const double pressure,
           const double fluidSoundSpeed) const {
}

}
}
