//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "GammaLawGas.cc"
#include "PhysicalConstants.hh"
#include "MKSUnits.hh"
#include "CGSUnits.hh"
#include "CosmologicalUnits.hh"
#include "SolarUnits.hh"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace Material {
    template class GammaLawGas<Dim<1>, PhysicalConstants<MKSUnits> >;
    template class GammaLawGas<Dim<2>, PhysicalConstants<MKSUnits> >;
    template class GammaLawGas<Dim<3>, PhysicalConstants<MKSUnits> >;
    template class GammaLawGas<Dim<1>, PhysicalConstants<CGSUnits> >;
    template class GammaLawGas<Dim<2>, PhysicalConstants<CGSUnits> >;
    template class GammaLawGas<Dim<3>, PhysicalConstants<CGSUnits> >;
    template class GammaLawGas<Dim<1>, PhysicalConstants<CosmologicalUnits> >;
    template class GammaLawGas<Dim<2>, PhysicalConstants<CosmologicalUnits> >;
    template class GammaLawGas<Dim<3>, PhysicalConstants<CosmologicalUnits> >;
    template class GammaLawGas<Dim<1>, PhysicalConstants<SolarUnits> >;
    template class GammaLawGas<Dim<2>, PhysicalConstants<SolarUnits> >;
    template class GammaLawGas<Dim<3>, PhysicalConstants<SolarUnits> >;
  }
}
