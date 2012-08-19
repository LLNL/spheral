//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "GruneisenEquationOfState.cc"
#include "Material/PhysicalConstants.hh"
#include "Material/MKSUnits.hh"
#include "Material/CGSUnits.hh"
#include "Material/SolarUnits.hh"

namespace Spheral {
  namespace SolidMaterial {
    template class GruneisenEquationOfState<Dim<1>, Material::PhysicalConstants<Material::MKSUnits> >;
    template class GruneisenEquationOfState<Dim<2>, Material::PhysicalConstants<Material::MKSUnits> >;
    template class GruneisenEquationOfState<Dim<3>, Material::PhysicalConstants<Material::MKSUnits> >;

    template class GruneisenEquationOfState<Dim<1>, Material::PhysicalConstants<Material::CGSUnits> >;
    template class GruneisenEquationOfState<Dim<2>, Material::PhysicalConstants<Material::CGSUnits> >;
    template class GruneisenEquationOfState<Dim<3>, Material::PhysicalConstants<Material::CGSUnits> >;

    template class GruneisenEquationOfState<Dim<1>, Material::PhysicalConstants<Material::SolarUnits> >;
    template class GruneisenEquationOfState<Dim<2>, Material::PhysicalConstants<Material::SolarUnits> >;
    template class GruneisenEquationOfState<Dim<3>, Material::PhysicalConstants<Material::SolarUnits> >;
  }
}
