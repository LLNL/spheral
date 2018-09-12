#include "PhysicsEvolvingMaterialLibrary.hh"
#include "Field/Field.hh"
#include "Utilities/bisectSearch.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/SpheralFunctions.hh"
#include "Utilities/DBC.hh"

#include <iostream>

namespace Spheral {

using std::min;
using std::max;
using std::abs;

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
PhysicsEvolvingMaterialLibrary<Dimension>::
PhysicsEvolvingMaterialLibrary(const double referenceDensity,
                               const double etamin,
                               const double etamax,
                               const PhysicalConstants& constants,
                               const double minimumPressure,
                               const double maximumPressure,
                               const MaterialPressureMinType minPressureType):
  Physics<Dimension>(),
  SolidEquationOfState<Dimension>(referenceDensity,
                                  etamin,
                                  etamax,
                                  constants,
                                  minimumPressure,
                                  maximumPressure,
                                  minPressureType),
  StrengthModel<Dimension>() {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
PhysicsEvolvingMaterialLibrary<Dimension>::
~PhysicsEvolvingMaterialLibrary() {
}

}
