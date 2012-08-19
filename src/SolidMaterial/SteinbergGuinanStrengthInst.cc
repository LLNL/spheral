//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SteinbergGuinanStrength.cc"
#include "Geometry/Dimension.hh"
#include "Material/PhysicalConstants.hh"
#include "Material/MKSUnits.hh"
#include "Material/CGSUnits.hh"
#include "Material/SolarUnits.hh"

using namespace Spheral::Material;

namespace Spheral {
  namespace SolidMaterial {
    template class SteinbergGuinanStrength<Dim<1>, PhysicalConstants<MKSUnits> >;
    template class SteinbergGuinanStrength<Dim<2>, PhysicalConstants<MKSUnits> >;
    template class SteinbergGuinanStrength<Dim<3>, PhysicalConstants<MKSUnits> >;

    template class SteinbergGuinanStrength<Dim<1>, PhysicalConstants<CGSUnits> >;
    template class SteinbergGuinanStrength<Dim<2>, PhysicalConstants<CGSUnits> >;
    template class SteinbergGuinanStrength<Dim<3>, PhysicalConstants<CGSUnits> >;

    template class SteinbergGuinanStrength<Dim<1>, PhysicalConstants<SolarUnits> >;
    template class SteinbergGuinanStrength<Dim<2>, PhysicalConstants<SolarUnits> >;
    template class SteinbergGuinanStrength<Dim<3>, PhysicalConstants<SolarUnits> >;
  }
}
