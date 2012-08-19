//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SteinbergGuinanLundStrength.cc"
#include "Geometry/Dimension.hh"
#include "Material/PhysicalConstants.hh"
#include "Material/MKSUnits.hh"
#include "Material/CGSUnits.hh"
#include "Material/SolarUnits.hh"

namespace Spheral {
  namespace SolidMaterial {
    template class SteinbergGuinanLundStrength<Dim<1>, PhysicalConstants<MKSUnits> >;
    template class SteinbergGuinanLundStrength<Dim<2>, PhysicalConstants<MKSUnits> >;
    template class SteinbergGuinanLundStrength<Dim<3>, PhysicalConstants<MKSUnits> >;

    template class SteinbergGuinanLundStrength<Dim<1>, PhysicalConstants<CGSUnits> >;
    template class SteinbergGuinanLundStrength<Dim<2>, PhysicalConstants<CGSUnits> >;
    template class SteinbergGuinanLundStrength<Dim<3>, PhysicalConstants<CGSUnits> >;

    template class SteinbergGuinanLundStrength<Dim<1>, PhysicalConstants<SolarUnits> >;
    template class SteinbergGuinanLundStrength<Dim<2>, PhysicalConstants<SolarUnits> >;
    template class SteinbergGuinanLundStrength<Dim<3>, PhysicalConstants<SolarUnits> >;
  }
}
