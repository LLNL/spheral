#include <iostream>
using namespace std;

#include "PhysicsEvolvingMaterialLibrary.hh"
#include "Field/Field.hh"
#include "Utilities/bisectSearch.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/SpheralFunctions.hh"
#include "Utilities/DBC.hh"

namespace Spheral {
namespace SolidMaterial {

using namespace std;
using std::min;
using std::max;
using std::abs;
using FieldSpace::Field;
using DataBaseSpace::DataBase;

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
PhysicsEvolvingMaterialLibrary<Dimension>::
PhysicsEvolvingMaterialLibrary(const Material::PhysicalConstants& constants,
                               const double minimumPressure,
                               const double maximumPressure,
                               const Material::MaterialPressureMinType minPressureType):
  Material::EquationOfState<Dimension>(constants,
                                       minimumPressure,
                                       maximumPressure,
                                       minPressureType),
  SolidMaterial::StrengthModel<Dimension>(),
  PhysicsSpace::Physics<Dimension>() {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
PhysicsEvolvingMaterialLibrary<Dimension>::
~PhysicsEvolvingMaterialLibrary() {
}

}
}

