//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "LinearPolynomialEquationOfState.cc"
#include "Geometry/Dimension.hh"
#include "Material/PhysicalConstants.hh"
#include "Material/MKSUnits.hh"
#include "Material/CGSUnits.hh"
#include "Material/SolarUnits.hh"

namespace Spheral {
  namespace SolidMaterial {
    template class LinearPolynomialEquationOfState<Dim<1>, Material::PhysicalConstants<Material::MKSUnits> >;
    template class LinearPolynomialEquationOfState<Dim<2>, Material::PhysicalConstants<Material::MKSUnits> >;
    template class LinearPolynomialEquationOfState<Dim<3>, Material::PhysicalConstants<Material::MKSUnits> >;

    template class LinearPolynomialEquationOfState<Dim<1>, Material::PhysicalConstants<Material::CGSUnits> >;
    template class LinearPolynomialEquationOfState<Dim<2>, Material::PhysicalConstants<Material::CGSUnits> >;
    template class LinearPolynomialEquationOfState<Dim<3>, Material::PhysicalConstants<Material::CGSUnits> >;

    template class LinearPolynomialEquationOfState<Dim<1>, Material::PhysicalConstants<Material::SolarUnits> >;
    template class LinearPolynomialEquationOfState<Dim<2>, Material::PhysicalConstants<Material::SolarUnits> >;
    template class LinearPolynomialEquationOfState<Dim<3>, Material::PhysicalConstants<Material::SolarUnits> >;
  }
}
