//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "MurnahanEquationOfState.cc"
#include "Geometry/Dimension.hh"
#include "Material/PhysicalConstants.hh"
#include "Material/MKSUnits.hh"
#include "Material/CGSUnits.hh"
#include "Material/SolarUnits.hh"

namespace Spheral {
  namespace SolidMaterial {
    template class MurnahanEquationOfState<Dim<1>, Material::PhysicalConstants<Material::MKSUnits> >;
    template class MurnahanEquationOfState<Dim<2>, Material::PhysicalConstants<Material::MKSUnits> >;
    template class MurnahanEquationOfState<Dim<3>, Material::PhysicalConstants<Material::MKSUnits> >;

    template class MurnahanEquationOfState<Dim<1>, Material::PhysicalConstants<Material::CGSUnits> >;
    template class MurnahanEquationOfState<Dim<2>, Material::PhysicalConstants<Material::CGSUnits> >;
    template class MurnahanEquationOfState<Dim<3>, Material::PhysicalConstants<Material::CGSUnits> >;

    template class MurnahanEquationOfState<Dim<1>, Material::PhysicalConstants<Material::SolarUnits> >;
    template class MurnahanEquationOfState<Dim<2>, Material::PhysicalConstants<Material::SolarUnits> >;
    template class MurnahanEquationOfState<Dim<3>, Material::PhysicalConstants<Material::SolarUnits> >;
  }
}
