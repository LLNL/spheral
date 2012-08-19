//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "NFWPotential.cc"
#include "Geometry/Dimension.hh"
#include "Material/PhysicalConstants.hh"
#include "Material/CGSUnits.hh"
#include "Material/MKSUnits.hh"
#include "Material/CosmologicalUnits.hh"

namespace Spheral {
  namespace PhysicsSpace {
    using namespace Material;
    template class NFWPotential<Dim<1>, PhysicalConstants<CGSUnits> >;
    template class NFWPotential<Dim<2>, PhysicalConstants<CGSUnits> >;
    template class NFWPotential<Dim<3>, PhysicalConstants<CGSUnits> >;

    template class NFWPotential<Dim<1>, PhysicalConstants<MKSUnits> >;
    template class NFWPotential<Dim<2>, PhysicalConstants<MKSUnits> >;
    template class NFWPotential<Dim<3>, PhysicalConstants<MKSUnits> >;
                                                      
    template class NFWPotential<Dim<1>, PhysicalConstants<CosmologicalUnits> >;
    template class NFWPotential<Dim<2>, PhysicalConstants<CosmologicalUnits> >;
    template class NFWPotential<Dim<3>, PhysicalConstants<CosmologicalUnits> >;
  }
}
